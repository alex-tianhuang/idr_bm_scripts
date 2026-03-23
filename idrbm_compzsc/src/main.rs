use anyhow::Error;
use bumpalo::{Bump, boxed::Box};
use clap::Parser;
use idrbm_core::box_dyn_in;
use std::{fs::File, io::Write, path::PathBuf};

use crate::utils::{read_variant_data, read_wt_data};
mod utils;
/// Program for computing some zscores by scaling wild-type data
/// using the mean and variance of some variant substitution data.
///
/// Currently extracts the last column of the given CSV.
#[derive(Parser)]
#[command(verbatim_doc_comment)]
struct Args {
    /// Input CSV of wild-type scores. Last column is where the score is assumed to be.
    wt_data: PathBuf,
    /// Input CSV of variant scores. Last column is where the score is assumed to be.
    variant_data: PathBuf,
    /// Output CSV file of regions and z-scores.
    output_file: PathBuf,
    /// When given, log errors from parsing the FASTA file
    /// or from invalid regions.
    ///
    /// Creates the file if does not exist,
    /// appends to the file if it does.
    #[arg(long)]
    log_errs_to: Option<PathBuf>,
}
fn main() -> Result<(), Error> {
    let Args {
        wt_data,
        variant_data,
        output_file,
        log_errs_to
    } = Args::try_parse()?;
    let arena = Bump::new();
    let mut err_out = alloc_error_handler(log_errs_to, &arena)?;
    let wt_data = read_wt_data(&wt_data, &arena, &mut *err_out)?;
    let variant_data = read_variant_data(&variant_data, &arena)?;

    let mut writer = csv::Writer::from_path(&output_file)?;
    writer.write_record(["ProteinID", "RegionID", "ZScore"])?;
    for (protein_id, protein_group) in variant_data {
        let raw = *wt_data.get(protein_id).ok_or_else(|| {
            Error::msg(format!(
                "could not find wild-type score for protein {}",
                protein_id
            ))
        })?;
        for &(region_id, [mean, std]) in protein_group {
            let z = (raw - mean) / std;
            writer.write_record([protein_id, region_id, &z.to_string()])?;
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