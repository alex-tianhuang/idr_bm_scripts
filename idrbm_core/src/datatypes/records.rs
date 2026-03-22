//! Module defining some csv record types.
//! 
//! See [`RegionCsvRecord`] and [`VariantCsvRecord`].

use crate::datatypes::aa_canonical_str;
/// 0-indexed bounds for a region.
/// 
/// Guaranteed that `start < stop`.
#[derive(Clone, Copy, Debug)]
pub struct RegionBounds {
    start: usize,
    stop: usize
}
/// The information for each row of a regions CSV.
pub struct RegionCsvRecord<'a> {
    pub protein_id: &'a str,
    pub region_id: &'a str,
    pub region: RegionBounds
}
/// The information for each row of a variant sequences CSV.
pub struct VariantCsvRecord<'a> {
    pub protein_id: &'a str,
    pub region_id: &'a str,
    pub variant_id: &'a str,
    pub variant_sequence: &'a aa_canonical_str
}
impl RegionBounds {
    /// Make a new [`RegionBounds`] if `start < stop`.
    pub fn new(start: usize, stop: usize) -> Option<Self> {
        if start >= stop {
            None
        } else {
            Some(Self { start, stop })
        }
    }
    /// The size of this region (`stop - start`).
    pub fn size(&self) -> usize {
        self.stop - self.start
    }
}