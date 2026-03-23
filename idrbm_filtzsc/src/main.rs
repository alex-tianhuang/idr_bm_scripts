use std::path::{Path, PathBuf};

use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use clap::Parser;
use csv::StringRecord;
use hashbrown::{DefaultHashBuilder, HashMap};
use idrbm_core::{datatypes::Grouped, utils::read_regions};
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
    let Args { regions, zscores, output_file, z } = Args::try_parse()?;
    if z <= 0.0 {
        return Err(Error::msg("z-score threshold must be positive"))
    }
    let arena = Bump::new();
    let regions = read_regions(&regions, &arena)?;
    let zscores = read_zscores(&zscores, &arena)?;

    let mut writer = csv::Writer::from_path(&output_file)?;
    writer.write_record(["ProteinID", "RegionID", "Start", "Stop"])?;

    for (protein_id, protein_group) in regions {
        let zscore_group = zscores.get(protein_id).ok_or_else(|| Error::msg(format!("failed to find zscore data for protein {}", protein_id)))?;
        for &(region_id, region) in protein_group {
            let &(_, zscore) = zscore_group.iter().find(|(rid, _)| *rid == region_id).ok_or_else(|| Error::msg(format!("failed to find zscore data for protein {} region {}", protein_id, region_id)))?;
            if zscore.abs() < z {
                continue;
            }
            writer.write_record([protein_id, region_id, &region.start().to_string(), &region.stop().to_string()])?;
        }
    }
    Ok(())
}
/// Will eventually refactor this to idrbm_core
/// but for now this reads zscores.
fn read_zscores<'a>(path: &Path, arena: &'a Bump) -> Result<Grouped<'a, &'a[(&'a str, f64)]>, Error> {
    let contents = std::fs::read(path)?;
    let mut reader = csv::Reader::from_reader(&*contents);
    let mut record = csv::StringRecord::new();
    let n_cols = reader.headers()?.len();
    if n_cols != 3 {
        return Err(Error::msg(format!("expected a `ProteinID`, `RegionID`, `ZScore` columns, got {} columns", n_cols)))
    }
    let local_arena = Bump::new();
    let mut local_order = Vec::new_in(&local_arena);
    let mut local_mapping: HashMap<
        &'a str,
        Vec<'_, (&'a str, f64)>,
        DefaultHashBuilder,
        &Bump,
    > = HashMap::new_in(&local_arena);
    while reader.read_record(&mut record)? {
        debug_assert_eq!(record.len(), n_cols);
        let [protein_id, region_id] = [0, 1].map(|i| record.get(i).unwrap());
        let zscore = record
            .get(2)
            .unwrap()
            .parse::<f64>()
            .map_err(|e| error(&e, &record, &contents))?;
        let protein_group = if let Some(group) = local_mapping.get_mut(protein_id) {
            group
        } else {
            let protein_id = &*arena.alloc_str(protein_id);
            local_order.push(protein_id);
            local_mapping
                .try_insert(protein_id, Vec::new_in(&local_arena))
                .unwrap()
        };
        let region_id = &*arena.alloc_str(region_id);
        protein_group.push((region_id, zscore))
    }
    let mut mapping = HashMap::with_capacity_in(local_mapping.len(), arena);
    for (&protein_id, protein_group) in local_mapping.iter() {
        // SAFETY: keys are from a map so they are unique
        unsafe {
            mapping.insert_unique_unchecked(
                protein_id,
                &*arena.alloc_slice_copy(&protein_group),
            )
        };
    }
    Ok(Grouped::new(arena.alloc_slice_copy(&local_order), mapping))
}
/// Add the line number and original line to the error
/// given by `reason`.
///
/// Helper for [`read_zscores`].
fn error(reason: &dyn std::fmt::Display, record: &StringRecord, data: &[u8]) -> Error {
    let position = record
        .position()
        .expect("position should be available if the record was successfully read");
    let data = &data[position.byte() as usize..];
    let next_newline = data.iter().position(|b| *b == b'\n').unwrap_or(data.len());
    let line = &data[..next_newline];
    Error::msg(format!(
        "unable to read line {} ({}): {}",
        position.line(),
        reason,
        String::from_utf8_lossy(line)
    ))
}