//! Module defining [`CompoundHeaderWriter`],
//! for taking a fixed number of IDs and joining them
//! using the '|' character.
/// The delimiter I use.
const DELIMITER: u8 = b'|';
/// Takes `N` IDs and joins them with '|' using `csv::Writer`.
/// 
/// Make one with [`CompoundHeaderWriter::new`].
pub struct CompoundHeaderWriter<const N: usize> {
    buf: Vec<u8>,
}
impl<const N: usize> CompoundHeaderWriter<N> {
    pub fn new() -> CompoundHeaderWriter<N> {
        CompoundHeaderWriter { buf: Vec::new() }
    }
    /// Get a header consisting of `parts` joined by the '|' character.
    pub fn construct_header(&mut self, parts: [&str; N]) -> &str {
        const EXPECT_MSG: &str = "failed to write to vec for some reason";
        let wtr = std::mem::take(&mut self.buf);
        let mut writer = csv::WriterBuilder::new().delimiter(DELIMITER).from_writer(wtr);
        writer.write_record(parts).expect(EXPECT_MSG);
        self.buf = writer.into_inner().expect(EXPECT_MSG);
        str::from_utf8(&self.buf).expect("unexpected non-UTF8 made by joining UTF-8 parts")
    }
}