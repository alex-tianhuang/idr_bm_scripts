//! File reading utilities.

use std::path::Path;

use anyhow::Error;
use bumpalo::Bump;
use csv::StringRecord;
use hashbrown::{DefaultHashBuilder, HashMap};
use idrbm_core::{
    csv::{CsvReader, DuplicateRule},
    datatypes::Grouped,
};

pub type WildTypeData<'a> = HashMap<&'a str, f64, DefaultHashBuilder, &'a Bump>;
/// Read the data from the wild-type data file at `path`.
pub fn read_wt_data<'a>(path: &Path, arena: &'a Bump) -> Result<WildTypeData<'a>, Error> {
    let contents = std::fs::read(path)?;
    CsvReader::new_with_at_least_n_cols::<2>(&|record, _| parse_f64_from_last_col(record))
        .collect_shallow_map(&contents, DuplicateRule::LastWins, arena)
}
pub type VariantData<'a> = Grouped<'a, &'a [(&'a str, [f64; 2])]>;
/// Read the means and standard deviations from the variant data file at `path`.
pub fn read_variant_data<'a>(path: &Path, arena: &'a Bump) -> Result<VariantData<'a>, Error> {
    let contents = std::fs::read(path)?;
    CsvReader::new_with_at_least_n_cols::<2>(&|record, _| parse_f64_from_last_col(record))
        .with_unpacker()
        .fold_regions_from_variants_with_compound_id(
            &contents,
            arena,
            || (0, [0.0, 0.0]),
            |(count, [sum, sum_sqr]), value| {
                *count += 1;
                *sum += value;
                *sum_sqr += value * value;
            },
            |(count, [sum, sum_sqr])| mean_and_sample_std(*count, *sum, *sum_sqr),
        )
}
/// Given a CSV row, get the last column and parse it as a number.
/// 
/// Helper for [`read_wt_data`] and [`read_wt_data`].
fn parse_f64_from_last_col(record: &StringRecord) -> Result<f64, Error> {
    let n_cols = record.len();
    Ok(record.get(n_cols - 1).unwrap().parse::<f64>()?)
}
/// Compute the mean and sample standard deviation.
///
/// Returns the `n-1` sample estimate.
pub fn mean_and_sample_std(count: usize, sum: f64, sum_of_squares: f64) -> [f64; 2] {
    let mean = sum / count as f64;
    let variance = (sum_of_squares - sum * mean) / (count - 1) as f64;
    [mean, variance.sqrt()]
}
