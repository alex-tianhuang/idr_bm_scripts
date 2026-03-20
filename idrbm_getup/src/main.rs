use crate::{checkpoint::Checkpoint, input::read_ids};
use anyhow::Error;
use bumpalo::Bump;
use clap::Parser;
use std::path::PathBuf;
mod checkpoint;
mod input;
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
    interval_retry: u64,
    /// When given, log errors from parsing Uniprot
    /// results to the given file instead of clouding
    /// `stderr`.
    ///
    /// Creates the file if does not exist,
    /// appends to the file if it does.
    #[arg(long)]
    log_parse_errs_to: Option<PathBuf>,
    /// When set, turn off the progress bar and spinner.
    #[arg(long, action)]
    disable_pbar: bool,
    /// When set, sets `--disable-pbar` and `--log-seq-errs-to`
    /// to the equivalent of `/dev/null`.
    #[arg(short, long, action)]
    quiet: bool,
}
/// See [`Args`] docs.
fn main() -> Result<(), Error> {
    let Args {
        input_file,
        output_file,
        interval_retry,
        log_parse_errs_to,
        disable_pbar,
        quiet,
    } = Args::try_parse()?;
    let arena = Bump::new();
    let raw_ids = read_ids(&input_file, &arena)?;
    let checkpoint = Checkpoint::load_or_empty(&output_file, &arena)?;
    let ids = checkpoint.retain(raw_ids);
    let enable_pbar = !(disable_pbar || quiet);
    let results_url = web::submit(ids, interval_retry, enable_pbar)?;
    smol::block_on(web::results(
        &output_file,
        &results_url,
        ids.len(),
        enable_pbar,
        &mut *web::alloc_error_handler(log_parse_errs_to, quiet, &arena)?,
    ))
}
