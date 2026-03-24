//! Module of csv reading tools.
//!
//! With little idea of how to organize this module
//! but a lot of CSV code with different deserialization styles,
//! [`CsvReader`] is the extremely general/difficult to document
//! utility and [`read_regions`] / [`read_variants`] are the ones
//! I will use 90% of the time.
mod compound_header;
mod duplicate_rule;
pub mod utils;
pub use crate::csv::compound_header::{CompoundHeaderReader, CompoundHeaderWriter};
use crate::{
    csv::utils::{collect_l2, collect_l3, finish_l2, finish_l3, fold_l2},
    datatypes::{RegionBounds, RegionMap, VariantMap, aa_canonical_str},
};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
pub use duplicate_rule::DuplicateRule;
use hashbrown::{DefaultHashBuilder, HashMap};
use std::{fmt::Debug, path::Path};
/// A sentinel value for [`CsvReader`]'s generic `N`.
const NO_UNPACKER: usize = 0;
/// An attempt at a flexible CSV deserializer that does not have
/// a signature that is too complicated (yes, I understand it still
/// looks too complicated).
///
/// Deserializes CSVs into [`Grouped`](crate::datatypes::Grouped) maps.
///
/// Usage
/// -----
/// Some barely organized usage tips.
///
/// Construct with [`CsvReader::new_with_at_least_n_cols`] to assert
/// that the CSV being read has a given number of columns.
///
/// Use [`CsvReader::with_unpacker`] to assert that the first column
/// will contain compound headers created by [`CompoundHeaderWriter`]
/// with `N_PARTS` parts to each header.
///
/// Use [`CsvReader::collect_regions`] or [`CsvReader::collect_variants`] to
/// gather a map of all region data/variant data, bucketed by protein ID (and region ID).
///
/// Some reasonable defaults using [`CsvReader::new_regions_reader`] and
/// [`CsvReader::new_variants_reader`].
pub struct CsvReader<'a, 'b, T, const N: usize = NO_UNPACKER> {
    deserializer: &'b dyn Fn(&csv::StringRecord, &'a Bump) -> Result<T, Error>,
    n_cols: usize,
    unpacker: Option<CompoundHeaderReader<N>>,
}
/// Read a regions CSV at the given `path`.
///
/// This file is expected have the following first four columns:
/// 1. ProteinID
/// 2. RegionID
/// 3. Start (int)
/// 4. Stop (int, greater than `Start`)
pub fn read_regions<'a>(
    path: &Path,
    duplicate_rule: DuplicateRule,
    arena: &'a Bump,
) -> Result<RegionMap<'a, RegionBounds>, Error> {
    let contents = std::fs::read(path)?;
    CsvReader::new_with_at_least_n_cols::<4>(&|record, _| {
        let [start, stop] = [2, 3].map(|i| {
            let field = record.get(i).unwrap();
            field.parse::<usize>()
        });
        let start = start?;
        let stop = stop?;
        RegionBounds::new(start, stop).ok_or_else(|| Error::msg("bad data, start >= stop"))
    })
    .collect_regions(&contents, duplicate_rule, arena)
}
/// Read a variant CSV at the given `path`.
///
/// This file is expected have the following first four columns:
/// 1. ProteinID
/// 2. RegionID
/// 3. VariantID
/// 4. Sequence
pub fn read_variants<'a>(
    path: &Path,
    duplicate_rule: DuplicateRule,
    arena: &'a Bump,
) -> Result<VariantMap<'a, &'a aa_canonical_str>, Error> {
    let contents = std::fs::read(path)?;
    CsvReader::new_with_at_least_n_cols::<4>(&|record, arena| {
        let s = aa_canonical_str::from_bytes(record.get(3).unwrap().as_bytes())?;
        Ok(aa_canonical_str::new(arena.alloc_slice_copy(s.as_slice())))
    })
    .collect_variants(&contents, duplicate_rule, arena)
}
impl<'a, 'b, T> CsvReader<'a, 'b, T> {
    /// Get a new CSV reader, asserting that there are at least
    /// `N_COLS` columns in the file that will be parsed.
    pub fn new_with_at_least_n_cols<const N_COLS: usize>(
        deserializer: &'b dyn Fn(&csv::StringRecord, &'a Bump) -> Result<T, Error>,
    ) -> Self {
        Self {
            deserializer,
            n_cols: N_COLS,
            unpacker: None,
        }
    }
    /// Parse a CSV with a compound header in the first column,
    /// as constructed by a [`CompoundHeaderWriter`] with `N_PARTS` parts.
    pub fn with_unpacker<const N_PARTS: usize>(self) -> CsvReader<'a, 'b, T, N_PARTS> {
        let Self {
            deserializer,
            n_cols,
            ..
        } = self;
        CsvReader {
            deserializer,
            n_cols,
            unpacker: Some(CompoundHeaderReader::new()),
        }
    }
    /// Collect and deserialize rows by
    /// the first column into an unordered hashmap.
    pub fn collect_shallow_map(
        self,
        contents: &[u8],
        duplicate_rule: DuplicateRule,
        arena: &'a Bump,
    ) -> Result<HashMap<&'a str, T, DefaultHashBuilder, &'a Bump>, Error> {
        let mut reader = csv::Reader::from_reader(contents);
        let mut record = csv::StringRecord::new();
        let Self {
            deserializer,
            n_cols,
            ..
        } = self;
        debug_assert!(
            n_cols >= 2,
            "should not read regions file with less than 2 expected columns"
        );
        check_at_least_n_cols(reader.headers()?.len(), n_cols)?;
        let mut mapping = HashMap::new_in(arena);
        while reader.read_record(&mut record)? {
            let protein_id = record.get(0).unwrap();
            let item = deserializer(&record, arena)
                .map_err(|reason| add_context(&reason, &record, contents))?;
            let k;
            if let Some((&protein_id, _)) = mapping.get_key_value(protein_id) {
                if matches!(duplicate_rule, DuplicateRule::FirstWins) {
                    continue;
                }
                k = protein_id;
            } else {
                k = &*arena.alloc_str(protein_id);
            }
            let _ = mapping.insert(k, item);
        }
        Ok(mapping)
    }
    /// Collect regions grouped by protein ID.
    ///
    /// This function will not check if regions are unique.
    pub fn collect_regions(
        self,
        contents: &[u8],
        duplicate_rule: DuplicateRule,
        arena: &'a Bump,
    ) -> Result<RegionMap<'a, T>, Error>
    where
        T: Copy + Debug,
    {
        let mut reader = csv::Reader::from_reader(contents);
        let mut record = csv::StringRecord::new();
        let Self {
            deserializer,
            n_cols,
            ..
        } = self;
        debug_assert!(
            n_cols >= 2,
            "should not read regions file with less than 2 expected columns"
        );
        check_at_least_n_cols(reader.headers()?.len(), n_cols)?;
        let local_arena = Bump::new();
        let mut local_mapping = HashMap::new_in(&local_arena);
        let mut local_order = Vec::new_in(&local_arena);
        while reader.read_record(&mut record)? {
            let [protein_id, region_id] = std::array::from_fn(|i| record.get(i).unwrap());
            let item = deserializer(&record, arena)
                .map_err(|reason| add_context(&reason, &record, contents))?;
            collect_l2(
                protein_id,
                region_id,
                item,
                duplicate_rule,
                &mut local_mapping,
                &mut local_order,
                &local_arena,
                arena,
            );
        }
        Ok(finish_l2(|item| *item, local_mapping, local_order, arena))
    }
    /// Collect variants grouped by protein ID and region ID.
    ///
    /// This function will not check if variants are unique.
    pub fn collect_variants(
        self,
        contents: &[u8],
        duplicate_rule: DuplicateRule,
        arena: &'a Bump,
    ) -> Result<VariantMap<'a, T>, Error>
    where
        T: Copy + Debug,
    {
        let mut reader = csv::Reader::from_reader(contents);
        let mut record = csv::StringRecord::new();
        let Self {
            deserializer,
            n_cols,
            ..
        } = self;
        debug_assert!(
            n_cols >= 3,
            "should not read regions file with less than 3 expected columns"
        );
        check_at_least_n_cols(reader.headers()?.len(), n_cols)?;
        let local_arena = Bump::new();
        let mut local_mapping = HashMap::new_in(&local_arena);
        let mut local_order = Vec::new_in(&local_arena);
        while reader.read_record(&mut record)? {
            let [protein_id, region_id, variant_id] =
                std::array::from_fn(|i| record.get(i).unwrap());
            let item = deserializer(&record, arena)
                .map_err(|reason| add_context(&reason, &record, contents))?;
            collect_l3(
                protein_id,
                region_id,
                variant_id,
                item,
                duplicate_rule,
                &mut local_mapping,
                &mut local_order,
                &local_arena,
                arena,
            );
        }
        Ok(finish_l3(|item| *item, local_mapping, local_order, arena))
    }
}
impl<'a, 'b, V> CsvReader<'a, 'b, V, 3> {
    // This is a spaghetti code tool used for one of my other crates
    // because it felt weird to implement it there.
    // Maybe this should have been done in a macro but that would
    // have taken even longer.
    // I digress, I'm not sure how to explain what this does
    // other than it's a fold operation over region-labelled data
    // that is not unique.
    pub fn fold_regions_from_variants_with_compound_id<Q, T>(
        self,
        contents: &[u8],
        arena: &'a Bump,
        default_items: impl Fn() -> Q,
        fold_items: impl Fn(&mut Q, V),
        finish_items: impl Fn(&Q) -> T,
    ) -> Result<RegionMap<'a, T>, Error>
    where
        Q: Copy + Debug,
    {
        let mut reader = csv::Reader::from_reader(contents);
        let mut record = csv::StringRecord::new();
        let Self {
            deserializer,
            n_cols,
            unpacker,
            ..
        } = self;
        let mut unpacker = unpacker
            .expect("cannot construct a CsvReader<..., N=3> without providing an unpacker!");
        debug_assert!(
            n_cols >= 2,
            "should not read regions file with less than 2 expected columns"
        );
        check_at_least_n_cols(reader.headers()?.len(), n_cols)?;
        let local_arena = Bump::new();
        let mut local_mapping = HashMap::new_in(&local_arena);
        let mut local_order = Vec::new_in(&local_arena);
        while reader.read_record(&mut record)? {
            let header = record.get(0).unwrap();
            let [protein_id, region_id, _] = unpacker
                .split_header(header)
                .map_err(|reason| add_context(&reason, &record, contents))?;
            let item = deserializer(&record, arena)
                .map_err(|reason| add_context(&reason, &record, contents))?;
            fold_l2(
                protein_id,
                region_id,
                item,
                &default_items,
                &fold_items,
                &mut local_mapping,
                &mut local_order,
                &local_arena,
                arena,
            );
        }
        Ok(finish_l2(finish_items, local_mapping, local_order, arena))
    }
}
impl<'a, 'b, T> CsvReader<'a, 'b, T, 3> {
    /// Like [`CsvReader::collect_variants`] but
    /// assumes a compound header in the first column.
    pub fn collect_variants_with_compound_id(
        self,
        contents: &[u8],
        duplicate_rule: DuplicateRule,
        arena: &'a Bump,
    ) -> Result<VariantMap<'a, T>, Error>
    where
        T: Copy + Debug,
    {
        let mut reader = csv::Reader::from_reader(contents);
        let mut record = csv::StringRecord::new();
        let Self {
            deserializer,
            n_cols,
            unpacker,
        } = self;
        let mut unpacker = unpacker
            .expect("cannot construct a CsvReader<..., N=3> without providing an unpacker!");
        debug_assert!(
            n_cols >= 1,
            "should not read compound-id variants file with less than 1 expected columns"
        );
        check_at_least_n_cols(reader.headers()?.len(), n_cols)?;
        let local_arena = Bump::new();
        let mut local_mapping = HashMap::new_in(&local_arena);
        let mut local_order = Vec::new_in(&local_arena);
        while reader.read_record(&mut record)? {
            let header = record.get(0).unwrap();
            let [protein_id, region_id, variant_id] = unpacker
                .split_header(header)
                .map_err(|reason| add_context(&reason, &record, contents))?;
            let item = deserializer(&record, arena)
                .map_err(|reason| add_context(&reason, &record, contents))?;
            collect_l3(
                protein_id,
                region_id,
                variant_id,
                item,
                duplicate_rule,
                &mut local_mapping,
                &mut local_order,
                &local_arena,
                arena,
            );
        }
        Ok(finish_l3(|item| *item, local_mapping, local_order, arena))
    }
}
impl<'a, 'b, T, const N: usize> CsvReader<'a, 'b, T, N> {}
/// Checks that the given number of columns is
/// at least `N` or fails otherwise.
fn check_at_least_n_cols(n_cols: usize, n_cols_expected: usize) -> Result<(), Error> {
    if n_cols >= n_cols_expected {
        Ok(())
    } else {
        Err(Error::msg(format!(
            "expected {} columns in CSV, found {}",
            n_cols_expected, n_cols
        )))
    }
}
/// Add the line number and original
/// line to the error given by `reason`.
fn add_context(
    reason: &dyn std::fmt::Display,
    record: &::csv::StringRecord,
    contents: &[u8],
) -> Error {
    let position = record
        .position()
        .expect("position should be available if the record was successfully read");
    let contents = &contents[position.byte() as usize..];
    let next_newline = contents
        .iter()
        .position(|b| *b == b'\n')
        .unwrap_or(contents.len());
    let line = &contents[..next_newline];
    Error::msg(format!(
        "unable to read line {} ({}): {}",
        position.line(),
        reason,
        String::from_utf8_lossy(line)
    ))
}
