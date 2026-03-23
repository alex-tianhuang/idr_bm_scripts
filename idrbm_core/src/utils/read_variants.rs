//! Module defining [`ParseVariants`], for parsing
//! a CSV of variants with four columns.

use std::{mem::ManuallyDrop, path::Path};

use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use csv::{Reader, StringRecord};
use hashbrown::HashMap;

use crate::datatypes::{VariantCsvRecord, Variants, aa_canonical_str};
/// Read a [`Variants`] from a variant sequences CSV at `path`.
///
/// These files contain the columns:
/// 1. ProteinID
/// 2. RegionID
/// 3. VariantID
/// 4. Variant sequence
///
/// All strings and arrays are allocated inside the given arena.
pub fn read_variants<'a>(path: &Path, arena: &'a Bump) -> Result<Variants<'a>, Error> {
    let contents = std::fs::read(path)?;
    let mut reader = ParseVariants::new(&contents);
    let mut record = StringRecord::new();

    let local_arena = Bump::new();
    let mut local_mapping = ManuallyDrop::new(HashMap::new_in(&local_arena));
    let mut local_order = Vec::new_in(&local_arena);
    while let Some(notification) = reader.next(&mut record) {
        let record = notification?;
        let protein_group = match local_mapping.get_mut(record.protein_id) {
            Some(group) => group,
            None => {
                let protein_id = &*arena.alloc_str(record.protein_id);
                local_order.push(protein_id);
                local_mapping
                    .try_insert(protein_id, Vec::new_in(&local_arena))
                    .unwrap()
            }
        };
        let (_, region_group) = match protein_group
            .iter_mut()
            .find(|(region_id, _)| *region_id == record.region_id)
        {
            Some(group) => group,
            None => {
                let region_id = &*arena.alloc_str(record.region_id);
                protein_group.push((region_id, Vec::new_in(&local_arena)));
                protein_group.last_mut().unwrap()
            }
        };
        let variant_id = &*arena.alloc_str(record.variant_id);
        let variant_sequence =
            aa_canonical_str::new(arena.alloc_slice_copy(record.variant_sequence.as_slice()));
        region_group.push((variant_id, variant_sequence));
    }
    let mut mapping = HashMap::with_capacity_in(local_mapping.len(), arena);
    for (protein_id, protein_group) in local_mapping.iter() {
        // SAFETY: keys are unique because they come from a map
        unsafe {
            mapping.insert_unique_unchecked(
                *protein_id,
                &*arena.alloc_slice_fill_iter(protein_group.iter().map(
                    |(region_id, region_group)| {
                        (*region_id, &*arena.alloc_slice_copy(&region_group))
                    },
                )),
            )
        };
    }
    Ok(Variants {
        mapping,
        order: arena.alloc_slice_copy(&local_order),
    })
}
/// Parser for variant sequence CSV files.
///
/// Helper for [`read_variants`].
struct ParseVariants<'a> {
    data: &'a [u8],
    inner: Reader<&'a [u8]>,
}
impl<'a> ParseVariants<'a> {
    /// Make a new [`ParseVariants`].
    pub fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            inner: Reader::from_reader(data),
        }
    }
    /// Read the next line into the given
    /// `buf` and attempt to parse it as a [`VariantCsvRecord`].
    ///
    /// Almost like an iterator but requires a `buf` to hold
    /// the result.
    ///
    /// Fails if the variant sequence is not canonical aminoacids.
    pub fn next<'b>(
        &mut self,
        buf: &'b mut StringRecord,
    ) -> Option<Result<VariantCsvRecord<'b>, Error>> {
        loop {
            match self.inner.read_record(buf) {
                Ok(true) => (),
                Ok(false) => return None,
                Err(e) => return Some(Err(Error::new(e))),
            }
            // empty stuff
            if buf.len() == 0 {
                continue;
            }
            return Some(self.parse_one_line(buf));
        }
    }
    /// Given a [`StringRecord`] that has been populated,
    /// attempt to parse a [`VariantCsvRecord`] from it.
    ///
    /// Fails if the variant sequence is not canonical aminoacids.
    ///
    /// Helper for [`ParseVariants::next`].
    fn parse_one_line<'b>(&self, record: &'b StringRecord) -> Result<VariantCsvRecord<'b>, Error> {
        if record.len() != 4 {
            return Err(error(
                &format_args!("expected four parts, got {}", record.len()),
                record,
                self.data,
            ));
        }
        let [protein_id, region_id, variant_id] = [0, 1, 2].map(|i| record.get(i).unwrap());
        let variant_sequence = aa_canonical_str::from_bytes(record.get(3).unwrap().as_bytes())
            .map_err(|e| error(&e, record, self.data))?;
        Ok(VariantCsvRecord {
            protein_id,
            region_id,
            variant_id,
            variant_sequence,
        })
    }
}

/// Add the line number and original line to the error
/// given by `reason`.
///
/// Helper for [`ParseVariants::parse_one_line`].
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
