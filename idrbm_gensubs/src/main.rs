use anyhow::Error;
use bumpalo::{Bump, boxed::Box};
use clap::Parser;
use idrbm_core::{
    box_dyn_in,
    csv::{CompoundHeaderWriter, DuplicateRule, read_regions, read_variants},
    utils::read_fasta_to_map,
};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};
/// Program for substituting variants into wild-type sequences.
#[derive(Parser)]
#[command(verbatim_doc_comment)]
struct Args {
    /// Input FASTA file of whole protein sequences.
    sequences: PathBuf,
    /// Input regions CSV file.
    regions: PathBuf,
    /// Input CSV file of variant sequences.
    variants: PathBuf,
    /// Output CSV of variant-substituted sequences.
    output_file: PathBuf,
    /// When given, log errors from parsing the FASTA file
    /// or from invalid regions.
    ///
    /// Creates the file if does not exist,
    /// appends to the file if it does.
    #[arg(long)]
    log_errs_to: Option<PathBuf>,
}
fn main() -> anyhow::Result<()> {
    let Args {
        sequences,
        regions,
        variants,
        output_file,
        log_errs_to,
    } = Args::try_parse()?;
    let arena = Bump::new();
    let mut err_out = alloc_error_handler(log_errs_to, &arena)?;
    let sequences = read_fasta_to_map(&sequences, &arena, &mut *err_out)?;
    let regions = read_regions(&regions, DuplicateRule::LastWins, &arena)?;
    let variants = read_variants(&variants, DuplicateRule::LastWins, &arena)?;

    let file = File::create(&output_file)?;
    let mut wtr = BufWriter::new(file);
    let mut header_writer = CompoundHeaderWriter::new();
    for (protein_id, protein_group) in variants {
        let whole_sequence = *sequences.get(protein_id).ok_or_else(|| {
            Error::msg(format!(
                "could not find whole protein sequence for protein {}",
                protein_id
            ))
        })?;
        let regions = regions.get(protein_id).ok_or_else(|| {
            Error::msg(format!(
                "could not find region bounds for protein {}",
                protein_id
            ))
        })?;
        for &(region_id, region_group) in protein_group {
            let (_, region) = regions
                .iter()
                .find(|rgn| rgn.0 == region_id)
                .ok_or_else(|| {
                    Error::msg(format!(
                        "could not find region bounds for protein {}, region {}",
                        protein_id, region_id
                    ))
                })?;
            if !region.is_in_bounds_of(whole_sequence) {
                writeln!(
                    err_out,
                    "region is out of bounds for protein {} (start={}, stop={}, seq.len()={}): {}",
                    region.start(),
                    region.stop(),
                    whole_sequence.len(),
                    protein_id,
                    region_id
                )?;
                continue;
            }
            let prefix = whole_sequence[..region.start()].as_str();
            let suffix = whole_sequence[region.stop()..].as_str();
            for &(variant_id, variant_sequence) in region_group {
                let compound_header =
                    header_writer.construct_header([protein_id, region_id, variant_id]);
                write!(
                    wtr,
                    ">{}\n{}{}{}\n",
                    compound_header,
                    prefix,
                    variant_sequence.as_str(),
                    suffix
                )?;
            }
        }
    }

    Ok(())
}

/// Helper function for callers of [`main`].
///
/// Construct a type erased thing that logs errors,
/// which can either be:
/// 1. `std::io::stderr`
/// 2. A `File` at a filepath given by the user.
fn alloc_error_handler(
    log_errs_to: Option<PathBuf>,
    arena: &Bump,
) -> Result<Box<'_, dyn Write>, Error> {
    match log_errs_to {
        Some(path) => {
            let writer = File::options().append(true).create(true).open(path)?;
            Ok(box_dyn_in!(Box::new_in(writer, arena) as Box<dyn Write>))
        }
        None => {
            let writer = std::io::stderr().lock();
            Ok(box_dyn_in!(Box::new_in(writer, arena) as Box<dyn Write>))
        }
    }
}
