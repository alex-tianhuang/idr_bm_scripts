//! Module defining [`RegionBounds`].

use crate::datatypes::aa_canonical_str;
/// 0-indexed bounds for a region.
/// 
/// Guaranteed that `start < stop`.
#[derive(Clone, Copy, Debug)]
pub struct RegionBounds {
    start: usize,
    stop: usize
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
    /// Utility method for checking that this regions
    /// is within bounds of the given `seq`.
    pub fn is_in_bounds_of(&self, seq: &aa_canonical_str) -> bool {
        (seq.len() >= self.stop) && (seq.len() >= self.start)
    }
    /// Inclusive lower bound.
    pub fn start(&self) -> usize {
        self.start
    }
    /// Exclusive upper bound.
    pub fn stop(&self) -> usize {
        self.stop
    }
}