use crate::datatypes::{Grouped, RegionBounds, RegionMap, VariantMap, aa_canonical_str};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use csv::{Reader, StringRecord};
use hashbrown::HashMap;
use std::{fmt::Debug, mem::ManuallyDrop, path::Path};

/// Read a regions CSV (containing `Start` and `Stop` columns as columns 3 and 4).
/// 
/// Otherwise see [`read_region_csv_template`].
pub fn read_regions<'a>(
    path: &Path,
    arena: &'a Bump,
) -> Result<RegionMap<'a, RegionBounds>, Error> {
    read_region_csv_template::<RegionBounds>(
        path,
        arena,
        |record| {
            let [start, stop] = [2, 3].map(|i| {
                let field = record.get(i).unwrap();
                field.parse::<usize>()
            });
            let start = start?;
            let stop = stop?;
            RegionBounds::new(start, stop).ok_or_else(|| Error::msg("bad data, start >= stop"))
        },
        assert_n_cols::<4>,
    )
}
/// Read a variant CSV (containing `Sequence` as column 4).
/// 
/// Otherwise see [`read_variant_csv_template`].
pub fn read_variants<'a>(
    path: &Path,
    arena: &'a Bump,
) -> Result<VariantMap<'a, &'a aa_canonical_str>, Error> {
    read_variant_csv_template::<&'a aa_canonical_str>(
        path,
        arena,
        |record| {
            let s = aa_canonical_str::from_bytes(record.get(3).unwrap().as_bytes())?;
            Ok(aa_canonical_str::new(arena.alloc_slice_copy(s.as_slice())))
        },
        assert_n_cols::<4>,
    )
}
/// Read region-labelled data from the given file.
///
/// These files contain the columns:
/// 1. ProteinID
/// 2. RegionID
/// and others.
///
/// The function signature is horrible
/// because it is extremely general.
///
/// Unconditionally panics if the number of columns is less
/// than 2 and `assert_n_cols` does not return an `Err`.
pub fn read_region_csv_template<'a, T: Debug + Copy>(
    path: &Path,
    arena: &'a Bump,
    deserializer: impl Fn(&StringRecord) -> Result<T, Error>,
    assert_n_cols: impl FnOnce(usize) -> Result<(), Error>,
) -> Result<RegionMap<'a, T>, Error> {
    let contents = std::fs::read(path)?;
    let mut reader = Reader::from_reader(&*contents);
    let n_cols = reader.headers()?.len();
    if let Err(e) = assert_n_cols(n_cols) {
        return Err(e);
    }
    if n_cols < 2 {
        panic!("cannot use read_region_csv_template with less than two columns")
    }
    let mut record = StringRecord::new();
    let local_arena = Bump::new();
    let mut local_mapping = ManuallyDrop::new(HashMap::new_in(&local_arena));
    let mut local_order = Vec::new_in(&local_arena);
    while reader.read_record(&mut record)? {
        let [protein_id, region_id] = std::array::from_fn(|i| {
            record
                .get(i)
                .expect("record length should be fixed to at least 2")
        });
        let item = deserializer(&record).map_err(|e| add_context(&e, &record, &contents))?;
        let items = match local_mapping.get_mut(protein_id) {
            Some(items) => items,
            None => {
                let protein_id = &*arena.alloc_str(protein_id);
                local_order.push(protein_id);
                let items = Vec::new_in(&local_arena);
                local_mapping.try_insert(protein_id, items).unwrap()
            }
        };
        let region_id = &*arena.alloc_str(region_id);
        items.push((region_id, item));
    }
    let mut mapping = HashMap::with_capacity_in(local_mapping.len(), arena);
    for (protein_id, items) in local_mapping.iter() {
        // SAFETY: keys are unique because they come from a map
        unsafe { mapping.insert_unique_unchecked(*protein_id, &*arena.alloc_slice_copy(&*items)) };
    }
    Ok(Grouped {
        mapping,
        order: arena.alloc_slice_copy(&local_order),
    })
}

/// Read variant-region-labelled data from the given `path`.
///  
/// These files contain the columns:
/// 1. ProteinID
/// 2. RegionID
/// 3. VariantID
/// and others.
///
/// The function signature is horrible
/// because it is extremely general.
///
/// Unconditionally panics if the number of columns is less
/// than 3 and `assert_n_cols` does not return an `Err`.
pub fn read_variant_csv_template<'a, T: Debug + Copy>(
    path: &Path,
    arena: &'a Bump,
    deserializer: impl Fn(&StringRecord) -> Result<T, Error>,
    assert_n_cols: impl FnOnce(usize) -> Result<(), Error>,
) -> Result<VariantMap<'a, T>, Error> {
    let contents = std::fs::read(path)?;
    let mut reader = Reader::from_reader(&*contents);
    let n_cols = reader.headers()?.len();
    if let Err(e) = assert_n_cols(n_cols) {
        return Err(e);
    }
    if n_cols < 2 {
        panic!("cannot use read_variant_csv_template with less than three columns")
    }
    let mut record = StringRecord::new();

    let local_arena = Bump::new();
    let mut local_mapping = ManuallyDrop::new(HashMap::new_in(&local_arena));
    let mut local_order = Vec::new_in(&local_arena);
    while reader.read_record(&mut record)? {
        let [protein_id, region_id, variant_id] = std::array::from_fn(|i| {
            record
                .get(i)
                .expect("record length should be fixed to at least 3")
        });
        let item = deserializer(&record).map_err(|e| add_context(&e, &record, &contents))?;
        let protein_group = match local_mapping.get_mut(protein_id) {
            Some(group) => group,
            None => {
                let protein_id = &*arena.alloc_str(protein_id);
                local_order.push(protein_id);
                local_mapping
                    .try_insert(protein_id, Vec::new_in(&local_arena))
                    .unwrap()
            }
        };
        let (_, items) = match protein_group.iter_mut().find(|(rid, _)| *rid == region_id) {
            Some(tuple) => tuple,
            None => {
                let region_id = &*arena.alloc_str(region_id);
                protein_group.push((region_id, Vec::new_in(&local_arena)));
                protein_group.last_mut().unwrap()
            }
        };
        let variant_id = &*arena.alloc_str(variant_id);
        items.push((variant_id, item));
    }
    let mut mapping = HashMap::with_capacity_in(local_mapping.len(), arena);
    for (protein_id, protein_group) in local_mapping.iter() {
        // SAFETY: keys are unique because they come from a map
        unsafe {
            mapping.insert_unique_unchecked(
                *protein_id,
                &*arena.alloc_slice_fill_iter(
                    protein_group
                        .iter()
                        .map(|(region_id, items)| (*region_id, &*arena.alloc_slice_copy(&*items))),
                ),
            )
        };
    }
    Ok(Grouped {
        mapping,
        order: arena.alloc_slice_copy(&local_order),
    })
}
/// Add the line number and original line to the error
/// given by `reason`.
fn add_context(reason: &dyn std::fmt::Display, record: &StringRecord, contents: &[u8]) -> Error {
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
/// Asserts that the given number of columns is `N`
/// or fails otherwise.
/// 
/// Helper for users of [`read_region_csv_template`]
/// and [`read_variant_csv_template`].
pub fn assert_n_cols<const N: usize>(n_cols: usize) -> Result<(), Error> {
    if n_cols == N {
        Ok(())
    } else {
        Err(Error::msg(format!(
            "expected {} columns in CSV, found {}",
            N, n_cols
        )))
    }
}
