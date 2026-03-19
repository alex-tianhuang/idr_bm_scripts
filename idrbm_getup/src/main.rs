use std::path::{Path, PathBuf};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use clap::Parser;
/// Program used to get sequences from Uniprot
/// given a file of Uniprot IDs. (hence GET UniProt, or GETUP).
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
    let input = read_input(&input_file, &arena)?;
    let checkpoint = read_checkpoint(&output_file, &arena)?;
    let ids = Vec::from_iter_in(input.iter().filter(|line| !checkpoint.contains(**line)).cloned(), &arena);
    getup(&ids, &output_file, max_connections)
}
fn read_input<'a>(path: &Path, arena: &'a Bump) -> Result<&'a [&'a str], Error> {
    todo!()
}
fn read_checkpoint<'a>(path: &Path, arena: &'a Bump) -> Result<std::mem::ManuallyDrop<hashbrown::HashSet<&'a str, hashbrown::DefaultHashBuilder, &'a Bump>>, Error> {
    todo!()
}
fn getup(ids: &[&str], output_file: &Path, max_connections: usize) -> Result<(), Error> {
    todo!()
}