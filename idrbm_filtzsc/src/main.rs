use anyhow::Error;
use bumpalo::Bump;
use clap::Parser;
use idrbm_core::{
    datatypes::RegionMap,
    utils::{assert_n_cols, read_region_csv_template, read_regions},
};
use std::path::{Path, PathBuf};
/// Program for computing some zscores by scaling wild-type data
/// using the mean and variance of some variant substitution data.
///
/// Currently extracts the last column of the given CSV.
#[derive(Parser)]
#[command(verbatim_doc_comment)]
struct Args {
    /// Input regions CSV.
    regions: PathBuf,
    /// Input CSV of region z-scores, made by `idrbm_compzsc`.
    zscores: PathBuf,
    /// Output CSV file of regions that pass a given threshold.
    output_file: PathBuf,
    /// Each region's absolute z-score must be at least this much
    /// to pass this filtering script.
    #[arg(short, default_value = "3")]
    z: f64,
}
fn main() -> Result<(), Error> {
    let Args {
        regions,
        zscores,
        output_file,
        z,
    } = Args::try_parse()?;
    if z <= 0.0 {
        return Err(Error::msg("z-score threshold must be positive"));
    }
    let arena = Bump::new();
    let regions = read_regions(&regions, &arena)?;
    let zscores = read_zscores(&zscores, &arena)?;

    let mut writer = csv::Writer::from_path(&output_file)?;
    writer.write_record(["ProteinID", "RegionID", "Start", "Stop"])?;

    for (protein_id, protein_group) in regions {
        let zscore_group = zscores.get(protein_id).ok_or_else(|| {
            Error::msg(format!(
                "failed to find zscore data for protein {}",
                protein_id
            ))
        })?;
        for &(region_id, region) in protein_group {
            let &(_, zscore) = zscore_group
                .iter()
                .find(|(rid, _)| *rid == region_id)
                .ok_or_else(|| {
                    Error::msg(format!(
                        "failed to find zscore data for protein {} region {}",
                        protein_id, region_id
                    ))
                })?;
            if zscore.abs() < z {
                continue;
            }
            writer.write_record([
                protein_id,
                region_id,
                &region.start().to_string(),
                &region.stop().to_string(),
            ])?;
        }
    }
    Ok(())
}
/// Load a file of region-labelled z-scores.
fn read_zscores<'a>(path: &Path, arena: &'a Bump) -> Result<RegionMap<'a, f64>, Error> {
    read_region_csv_template(
        path,
        arena,
        |record| {
            let value = record.get(2).unwrap().parse::<f64>()?;
            Ok(value)
        },
        assert_n_cols::<3>,
    )
}
