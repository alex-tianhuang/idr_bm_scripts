//! File reading utilities.

use std::{io::Write, path::Path};

use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use csv::{Reader, StringRecord};
use hashbrown::{DefaultHashBuilder, HashMap};
use idrbm_core::{datatypes::Grouped, utils::CompoundHeaderReader};

pub type WildTypeData<'a> = HashMap<&'a str, f64, DefaultHashBuilder, &'a Bump>;
/// Read the data from the wild-type data file at `path`.
pub fn read_wt_data<'a>(
    path: &Path,
    arena: &'a Bump,
    err_out: &mut dyn Write,
) -> Result<WildTypeData<'a>, Error> {
    let contents = std::fs::read(path)?;
    let mut reader = Reader::from_reader(&*contents);
    let mut record = StringRecord::new();
    let n_cols = reader.headers()?.len();
    if n_cols < 2 {
        return Err(Error::msg("expected at least two columns"));
    }
    let mut data = HashMap::new_in(arena);
    while reader.read_record(&mut record)? {
        debug_assert_eq!(record.len(), n_cols);
        let protein_id = record.get(0).unwrap();
        let value = record
            .get(n_cols - 1)
            .unwrap()
            .parse::<f64>()
            .map_err(|e| error(&e, &record, &contents))?;
        let k;
        if let Some((&protein_id, _)) = data.get_key_value(protein_id) {
            writeln!(err_out, "duplicate entry, evicting oldest: {}", protein_id)?;
            k = protein_id;
        } else {
            k = &*arena.alloc_str(protein_id);
        }
        let _ = data.insert(k, value);
    }
    Ok(data)
}
pub type VariantData<'a> = Grouped<'a, &'a [(&'a str, [f64; 2])]>;
/// Read the means and standard deviations from the variant data file at `path`.
pub fn read_variant_data<'a>(
    path: &Path,
    arena: &'a Bump,
) -> Result<VariantData<'a>, Error> {
    let contents = std::fs::read(path)?;
    let mut reader = Reader::from_reader(&*contents);
    let mut header_reader = <CompoundHeaderReader<3>>::new();
    let mut record = StringRecord::new();
    let n_cols = reader.headers()?.len();
    if n_cols < 2 {
        return Err(Error::msg("expected at least two columns"));
    }
    let local_arena = Bump::new();
    let mut local_order = Vec::new_in(&local_arena);
    let mut local_mapping: HashMap<
        &'a str,
        Vec<'_, (&'a str, (usize, f64, f64))>,
        DefaultHashBuilder,
        &Bump,
    > = HashMap::new_in(&local_arena);
    while reader.read_record(&mut record)? {
        debug_assert_eq!(record.len(), n_cols);
        let compound_id = record.get(0).unwrap();
        let [protein_id, region_id] = match header_reader.split_header(compound_id) {
            Ok([protein_id, region_id, _]) => [protein_id, region_id],
            Err(e) => return Err(error(&e, &record, &contents)),
        };
        let value = record
            .get(n_cols - 1)
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
        match protein_group.iter_mut().find(|(rid, _)| *rid == region_id) {
            Some((_, (count, sum, sum_sqr))) => {
                *count += 1;
                *sum += value;
                *sum_sqr += value * value;
            }
            None => {
                let region_id = &*arena.alloc_str(region_id);
                protein_group.push((region_id, (1, value, value * value)))
            }
        }
    }
    let mut mapping = HashMap::with_capacity_in(local_mapping.len(), arena);
    for (&protein_id, protein_group) in local_mapping.iter() {
        // SAFETY: keys are from a map so they are unique
        unsafe {
            mapping.insert_unique_unchecked(
                protein_id,
                &*arena.alloc_slice_fill_iter(protein_group.iter().map(
                    |&(region_id, (count, sum, sum_sqr))| {
                        (region_id, mean_and_sample_std(count, sum, sum_sqr))
                    },
                )),
            )
        };
    }
    Ok(VariantData::new(
        arena.alloc_slice_copy(&local_order),
        mapping,
    ))
}
/// Add the line number and original line to the error
/// given by `reason`.
///
/// Helper for [`read_wt_data`] and [`read_variant_data`].
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
/// Compute the mean and sample standard deviation.
///
/// Returns the `n-1` sample estimate.
pub fn mean_and_sample_std(count: usize, sum: f64, sum_of_squares: f64) -> [f64; 2] {
    let mean = sum / count as f64;
    let variance = (sum_of_squares - sum * mean) / (count - 1) as f64;
    [mean, variance.sqrt()]
}
