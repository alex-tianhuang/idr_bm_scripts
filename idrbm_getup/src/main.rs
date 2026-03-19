use std::path::{Path, PathBuf};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use clap::Parser;
use crate::{checkpoint::Checkpoint, input::read_ids};
mod input;
mod checkpoint;

/// A program used to get sequences from the Uniprot
/// webservice online given a file of Uniprot IDs.
#[derive(Parser)]
#[command(verbatim_doc_comment)]
pub struct Args {
    /// Input text file of Uniprot IDs, one per line.
    input_file: PathBuf,
    /// Output FASTA file of canonical Uniprot sequences.
    output_file: PathBuf,
    /// Maximum number of concurrent requests to make.
    #[arg(short = 'C', long, default_value = "50")]
    max_connections: usize
}
fn main() -> Result<(), Error> {
    let Args { input_file, output_file, max_connections } = Args::parse();
    let arena = Bump::new();
    let raw_ids = read_ids(&input_file, &arena)?;
    let checkpoint = Checkpoint::load_or_empty(&output_file, &arena)?;
    let ids = checkpoint.retain(raw_ids);
    getup(&ids, &output_file, max_connections)
}
fn read_checkpoint<'a>(path: &Path, arena: &'a Bump) -> Result<std::mem::ManuallyDrop<hashbrown::HashSet<&'a str, hashbrown::DefaultHashBuilder, &'a Bump>>, Error> {
    todo!()
}
fn getup(ids: &[&str], output_file: &Path, max_connections: usize) -> Result<(), Error> {
    todo!()
}