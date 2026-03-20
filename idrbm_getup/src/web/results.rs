//! Module defining [`results`], a function which
//! downloads results from `https://rest.uniprot.org/idmapping/results`.
use anyhow::Error;
use bumpalo::{Bump, boxed::Box};
use idrbm_core::{box_dyn_in, datatypes::aa_canonical_str};
use indicatif::{ProgressBar, ProgressStyle};
use reqwest::{Client, Response};
use serde::Deserialize;
use serde_json::value::RawValue;
use tokio::{io::AsyncWriteExt, spawn, task::yield_now};
use std::{fs::File, io::Write, path::{Path, PathBuf}};

/// Helper function for callers of [`results`].
///
/// Construct a type erased thing that logs errors,
/// which can either be:
/// 1. `std::io::stderr`
/// 2. A `File` at a filepath given by the user.
/// 3. Nothing on `--quiet` mode.
pub fn alloc_error_handler(
    log_parse_errs_to: Option<PathBuf>,
    quiet_override: bool,
    arena: &Bump,
) -> Result<Box<'_, dyn Write>, Error> {
    if quiet_override {
        let writer = std::io::sink();
        return Ok(box_dyn_in!(Box::new_in(writer, arena) as Box<dyn Write>))
    }
    match log_parse_errs_to {
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

/// Fetch Uniprot idmapping results from the given `results_url`.
/// 
/// `results_url` is the URL returned by [`super::submit::submit`].
pub async fn results(
    output_path: &Path,
    results_url: &str,
    num_seqs_expected: usize,
    enable_pbar: bool,
    err_out: &mut dyn Write,
) -> anyhow::Result<()> {
    let mut file = tokio::fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(output_path)
        .await?;
    let client = Client::new();
    let mut request = Some(spawn(client.get(results_url).send()));
    let pbar = enable_pbar.then(|| pbar(num_seqs_expected));
    while let Some(req) = request {
        let response = req.await??.error_for_status()?;
        let next_page_url = extract_next_page_url(&response);
        request = next_page_url.map(|url| spawn(client.get(url).send()));
        let bytes = response.bytes().await?;
        if let Ok(ResultsResponse { results }) = serde_json::from_slice(&bytes) {
            // hopefully give `request` a chance to start before doing a bunch of deserialization
            yield_now().await;
            let results = <Vec<&RawValue>>::deserialize(results)?;
            write_results(&mut file, &results, err_out).await?;
            pbar.as_ref().map(|pbar| pbar.inc(results.len() as u64));
        }
    }
    pbar.as_ref().map(ProgressBar::abandon);
    Ok(())
}
/// Get the URL with the next page of results.
fn extract_next_page_url(response: &Response) -> Option<&str> {
    response.headers().get(reqwest::header::LINK).and_then(|h| {
        let header = h.to_str().ok()?;
        header.split(',').find_map(|part| {
            let part = part.trim();
            if part.ends_with("rel=\"next\"") {
                let start = part.find('<')? + 1;
                let end = part.find('>')?;
                Some(&part[start..end])
            } else {
                None
            }
        })
    })
}
/// Partial response @ `https://rest.uniprot.org/idmapping/results`.
/// 
/// See also [`write_results`] and [`extract_id_and_sequence`].
#[derive(Deserialize)]
struct ResultsResponse<'a> {
    #[serde(borrow)]
    // Vec<&'a RawValue>
    results: &'a RawValue,
}
/// Write the given results to `file`.
/// 
/// This function fails if `err_out` cannot be written to
/// or if `file` cannot be written to.
/// 
/// The function will not fail due to any of the `RawValue`s
/// failing to deserialize.
async fn write_results(
    file: &mut tokio::fs::File,
    results: &[&RawValue],
    err_out: &mut dyn Write,
) -> anyhow::Result<()> {
    let mut buf = Vec::new();
    for item in results.iter() {
        let parsed;
        match extract_id_and_sequence(item) {
            Ok((uniprot_id, sequence)) => {
                if let Err(e) = aa_canonical_str::from_bytes(sequence.as_bytes()) {
                    writeln!(
                        err_out,
                        "failed to process Uniprot entry ({}): {}",
                        e, uniprot_id
                    )?;
                    continue;
                }
                if sequence.is_empty() {
                    writeln!(
                        err_out,
                        "failed to process Uniprot entry (empty sequence): {}",
                        uniprot_id
                    )?;
                    continue;
                }
                parsed = (uniprot_id, sequence);
            }
            Err(e) => {
                writeln!(err_out, "{}", e)?;
                continue;
            }
        }
        let (uniprot_id, sequence) = parsed;
        write!(&mut buf, ">{}\n{}\n", uniprot_id, sequence)?;
    }
    file.write_all(&buf).await?;
    Ok(())
}
/// Deserialize a [`RawValue`] into a `(uniprot_id, sequence)` pair.
/// 
/// Attempts to access the `from` field on failure for a better error message.
fn extract_id_and_sequence(item: &RawValue) -> anyhow::Result<(&str, &str)> {
    #[derive(Deserialize)]
    struct SequenceEntry<'a> {
        value: &'a str,
    }
    #[derive(Deserialize)]
    struct ToEntry<'a> {
        #[serde(borrow)]
        sequence: SequenceEntry<'a>,
    }
    #[derive(Deserialize)]
    struct ResultEntry<'a> {
        from: &'a str,
        to: ToEntry<'a>,
    }
    #[derive(Deserialize)]
    struct FromEntry<'a> {
        from: &'a str,
    }
    let entry =
        <ResultEntry>::deserialize(item).map_err(|e| match <FromEntry>::deserialize(item) {
            Ok(entry) => Error::msg(format!(
                "failed to process Uniprot entry ({}): {}",
                e, entry.from
            )),
            Err(e) => e.into(),
        })?;
    Ok((entry.from, entry.to.sequence.value))
}
/// Progress bar for the [`results`] function.
fn pbar(n: usize) -> ProgressBar {
    ProgressBar::new(n as u64)
        .with_style(
            ProgressStyle::default_bar()
                .template("{msg} [{bar:40}] {pos}/{len}, ETA {eta}")
                .unwrap(),
        )
        .with_message("downloading sequences from Uniprot")
}
