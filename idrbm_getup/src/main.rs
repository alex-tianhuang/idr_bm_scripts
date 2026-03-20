use std::path::{Path, PathBuf};
use anyhow::Error;
use bumpalo::{Bump};
use clap::Parser;
use crate::{checkpoint::Checkpoint, input::read_ids};
mod input;
mod checkpoint;
mod web;

/// A program used to get sequences from the Uniprot
/// webservice online given a file of Uniprot IDs.
#[derive(Parser)]
#[command(verbatim_doc_comment)]
pub struct Args {
    /// Input text file of Uniprot IDs, one per line.
    input_file: PathBuf,
    /// Output FASTA file of canonical Uniprot sequences.
    output_file: PathBuf,
    /// Interval (in seconds) between polls to Uniprot
    /// to check whether your job is done.
    #[arg(short = 'I', long, default_value = "1")]
    interval_retry: u64
}
fn main() -> Result<(), Error> {
    let Args { input_file, output_file, interval_retry } = Args::try_parse()?;
    let arena = Bump::new();
    let raw_ids = read_ids(&input_file, &arena)?;
    let checkpoint = Checkpoint::load_or_empty(&output_file, &arena)?;
    let ids = checkpoint.retain(raw_ids);
    web::submit(ids, interval_retry, true)?;
    todo!()
}
fn getup(ids: &[&str], output_file: &Path, max_connections: usize) -> Result<(), Error> {
    todo!()
}