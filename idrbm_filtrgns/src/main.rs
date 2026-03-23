use std::{fs::File, io::{BufWriter, Write}, path::PathBuf};

use bumpalo::{Bump, boxed::Box};
use clap::Parser;
use anyhow::Error;
use idrbm_core::{box_dyn_in, utils::{read_fasta, read_regions}};
/// Program for filtering a FASTA file
/// for sequences with regions in the given regions file.
#[derive(Parser)]
#[command(verbatim_doc_comment)]
struct Args {
    /// Input FASTA.
    sequences: PathBuf,
    /// Input regions CSV.
    regions: PathBuf,
    /// Output FASTA with only the protein IDs in the given regions file.
    output_file: PathBuf,
    /// When given, log errors from parsing the FASTA file.
    ///
    /// Creates the file if does not exist,
    /// appends to the file if it does.
    #[arg(long)]
    log_errs_to: Option<PathBuf>,
}
fn main() -> Result<(), Error>{
    let Args { sequences, regions, output_file, log_errs_to } = Args::try_parse()?;
    let arena = Bump::new();
    let sequences = read_fasta(&sequences, &arena, &mut * alloc_error_handler(log_errs_to, &arena)?)?;
    let regions = read_regions(&regions, &arena)?;
    let file = File::create(&output_file)?;
    let mut writer = BufWriter::new(file);
    for entry in sequences {
        if regions.contains_key(entry.header) {
            write!(writer, ">{}\n{}\n", entry.header, entry.sequence.as_str())?;
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