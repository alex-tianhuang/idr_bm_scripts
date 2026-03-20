use anyhow::Error;
use indicatif::{ProgressBar, ProgressStyle};
use reqwest::blocking::Client;
use serde::Deserialize;
use std::{io::Read, time::Duration};

/// Submit and poll a Uniprot job to completion,
/// without downloading the results.
/// 
/// Returns the job ID of the completed job.
pub fn submit(ids: &[&str], interval_retry: u64, enable_spinner: bool) -> anyhow::Result<String> {
    let client = Client::new();
    let bytes = post_run(&client, ids)?;
    let RunResponse { job_id } = serde_json::from_slice(&bytes)?;
    let spinner = enable_spinner.then(pbar_spinner);
    poll_status_until_done(&client, job_id, interval_retry, spinner)?;
    let bytes = get_details(&client, job_id)?;
    let DetailsResponse { redirect_url } = serde_json::from_slice(&bytes)?;
    Ok(redirect_url.to_owned())
}
/// The URL for Uniprot's idmapping service.
const API_URL: &str = "https://rest.uniprot.org/idmapping";
/// Response @ `https://rest.uniprot.org/idmapping/run`
#[derive(Deserialize)]
struct RunResponse<'a> {
    #[serde(rename = "jobId")]
    job_id: &'a str,
}
/// Post to `https://rest.uniprot.org/idmapping/run`
/// and return the response as bytes.
fn post_run<'a>(client: &Client, ids: &[&str]) -> anyhow::Result<Vec<u8>> {
    let mut response = client
        .post(format!("{}/run", API_URL))
        .form(&[
            ("from", "UniProtKB_AC-ID"),
            ("to", "UniProtKB"),
            ("ids", &ids.join(",")),
        ])
        .send()?
        .error_for_status()?;
    let mut bytes = Vec::new();
    response.read_to_end(&mut bytes)?;
    Ok(bytes)
}
/// Response @ `https://rest.uniprot.org/idmapping/status`
#[derive(Deserialize)]
struct StatusResponse<'a> {
    #[serde(rename = "jobStatus")]
    job_status: Option<&'a str>,
}
/// Repeatedly poll `https://rest.uniprot.org/idmapping/status`
/// until it signals that the underlying job is done.
fn poll_status_until_done(
    client: &Client,
    job_id: &str,
    interval_retry: u64,
    spinner: Option<ProgressBar>
) -> anyhow::Result<()> {
    loop {
        let bytes = poll_status(&client, job_id)?;
        let StatusResponse { job_status } = serde_json::from_slice(&bytes)?;

        match job_status.as_deref() {
            Some("RUNNING") | Some("NEW") => {
                std::thread::sleep(Duration::from_secs(interval_retry));
                spinner.as_ref().map(|spinner| spinner.inc(1));
                continue;
            }
            Some("FINISHED") | Some("SUCCESS") | None => return Ok(()),
            Some(other) => {
                return Err(Error::msg(format!(
                    "job ended unexpectedly with status: {}",
                    other
                )));
            }
        }
    }
}
/// Send a get req to `https://rest.uniprot.org/idmapping/status`
/// and return the response as bytes.
fn poll_status<'a>(client: &Client, job_id: &str) -> anyhow::Result<Vec<u8>> {
    let mut response = client
        .get(format!("{}/status/{}", API_URL, job_id))
        .send()?
        .error_for_status()?;
    let mut bytes = Vec::new();
    response.read_to_end(&mut bytes)?;
    Ok(bytes)
}
/// Response @ `https://rest.uniprot.org/idmapping/details`
#[derive(Deserialize)]
struct DetailsResponse<'a> {
    #[serde(rename = "redirectURL")]
    redirect_url: &'a str,
}
/// Send a get req to `https://rest.uniprot.org/idmapping/details`
/// and return the response as bytes.
fn get_details(client: &Client, job_id: &str) -> anyhow::Result<Vec<u8>> {
    let mut response = client
        .get(format!("{}/idmapping/details/{}", API_URL, job_id))
        .send()?
        .error_for_status()?;
    let mut bytes = Vec::new();
    response.read_to_end(&mut bytes)?;
    Ok(bytes)
}
/// Spinner for [`submit`].
fn pbar_spinner() -> ProgressBar {
    ProgressBar::new_spinner().with_style(ProgressStyle::default_spinner()
        .template("{msg} [# requests: {pos}] {spinner}")
        .unwrap()).with_position(1).with_message("polling <uniprot>/idmapping/status")
}