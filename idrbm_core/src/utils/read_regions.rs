//! Module defining [`ParseRegions`] and [`read_regions`],
//! for parsing a CSV of regions with four columns.

use std::{mem::ManuallyDrop, path::Path};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use csv::{Reader, StringRecord};
use hashbrown::HashMap;

use crate::datatypes::{RegionBounds, RegionCsvRecord, Regions};
/// Read a [`Regions`] from a regions CSV at `path`.
/// 
/// These files contain the columns:
/// 1. ProteinID
/// 2. RegionID
/// 3. Start (number)
/// 4. Stop (number)
/// 
/// All strings and arrays are allocated inside the given arena.
pub fn read_regions<'a>(path: &Path, arena: &'a Bump) -> Result<Regions<'a>, Error> {
    let contents = std::fs::read(path)?;
    let mut reader = ParseRegions::new(&contents);
    let mut record = StringRecord::new();

    let local_arena = Bump::new();
    let mut local_mapping = ManuallyDrop::new(HashMap::new_in(&local_arena));
    let mut local_order = Vec::new_in(&local_arena);
    while let Some(notification) = reader.next(&mut record) {
        let record = notification?;
        let group = match local_mapping.get_mut(record.protein_id) {
            Some(group) => {
                group
            },
            None => {
                let protein_id = &*arena.alloc_str(record.protein_id);
                local_order.push(protein_id);
                local_mapping.try_insert(protein_id, Vec::new_in(&local_arena)).unwrap()
            }
        };
        let region_id = &*arena.alloc_str(record.region_id);
        group.push((region_id, record.region))
    }
    let mut mapping = HashMap::with_capacity_in(local_mapping.len(), arena);
    for (protein_id, group) in local_mapping.iter() {
        // SAFETY: keys are unique because they come from a map
        unsafe { mapping.insert_unique_unchecked(*protein_id, &*arena.alloc_slice_copy(&group)) };
    }
    Ok(Regions {
        mapping,
        order: arena.alloc_slice_copy(&local_order)
    })
}
/// Parser for regions CSV files.
/// 
/// Helper struct for [`read_regions`].
struct ParseRegions<'a> {
    data: &'a [u8],
    inner: Reader<&'a [u8]>,
}
impl<'a> ParseRegions<'a> {
    /// Make a new [`ParseRegions`].
    pub fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            inner: Reader::from_reader(data),
        }
    }
    /// Read the next line into the given
    /// `buf` and attempt to parse it as a [`RegionCsvRecord`].
    /// 
    /// Almost like an iterator but requires a `buf` to hold
    /// the result.
    pub fn next<'b>(
        &mut self,
        buf: &'b mut StringRecord,
    ) -> Option<Result<RegionCsvRecord<'b>, Error>> {
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
            return Some(self.parse_one_line(buf))
        }
    }
    /// Given a [`StringRecord`] that has been populated,
    /// attempt to parse a [`RegionCsvRecord`] from it.
    /// 
    /// Helper for [`ParseRegions::next`].
    fn parse_one_line<'b>(&self, record: &'b StringRecord) -> Result<RegionCsvRecord<'b>, Error> {
        if record.len() != 4 {
            return Err(error(
                &format_args!("expected four parts, got {}", record.len()),
                record,
                self.data,
            ));
        }
        let [protein_id, region_id] = [0, 1].map(|i| record.get(i).unwrap());
        let [start, stop] = [2, 3].map(|i| {
            let field = record.get(i).unwrap();
            field
                .parse::<usize>()
                .map_err(|e| error(&e, record, self.data))
        });
        let start = start?;
        let stop = stop?;
        let region = RegionBounds::new(start, stop).ok_or_else(|| {
            error(&"bad data, start >= stop", record, self.data)
        })?;
        Ok(RegionCsvRecord { protein_id, region_id, region })
    }
}

/// Add the line number and original line to the error
/// given by `reason`.
/// 
/// Helper for [`ParseRegions::parse_one_line`].
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
