/// The information for each row of a regions CSV.
/// 
/// A type guarantee is that `start` < `stop`.
pub struct RegionCsvRecord<'a> {
    pub protein_id: &'a str,
    pub region_id: &'a str,
    start: usize,
    stop: usize
}
impl<'a> RegionCsvRecord<'a> {
    /// Make a new [`RegionCsvRecord`] if `start < stop`.
    pub fn new(protein_id: &'a str, region_id: &'a str, start: usize, stop: usize) -> Option<Self> {
        if start >= stop {
            None
        } else {
            Some(Self { protein_id, region_id, start, stop })
        }
    }
    /// Inclusive start coordinate.
    pub fn start(&self) -> usize {
        self.start
    }
    /// Exclusive stop coordinate.
    pub fn stop(&self) -> usize {
        self.stop
    }
    /// Length of the region.
    /// `self.stop - self.start`.
    pub fn region_len(&self) -> usize {
        self.stop - self.start
    }
}