use anyhow::Error;
use bumpalo::Bump;
use clap::Parser;
use std::path::PathBuf;

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
}
fn main() -> Result<(), Error> {
    let Args {
        wt_data,
        variant_data,
        output_file,
    } = Args::try_parse()?;
    let arena = Bump::new();
    let wt_data = read_wt_data(&wt_data, &arena)?;
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
