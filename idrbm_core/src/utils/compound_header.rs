//! Module for taking a fixed number of IDs
//! and joining them using the '|' character,
//! or doing the opposite (decoding the IDs
//! from a compound ID).

use anyhow::Error;
use csv::StringRecord;
/// The delimiter I use.
const DELIMITER: u8 = b'|';
/// Takes `N` IDs and joins them with '|' using `csv::Writer`.
///
/// Make one with [`CompoundHeaderWriter::new`].
pub struct CompoundHeaderWriter<const N: usize> {
    buf: Vec<u8>,
}
/// Takes an ID made from a [`CompoundHeaderWriter<N>`]
/// and returns the parts that made it.
pub struct CompoundHeaderReader<const N: usize> {
    buf: StringRecord,
}
impl<const N: usize> CompoundHeaderWriter<N> {
    pub fn new() -> Self {
        Self { buf: Vec::new() }
    }
    /// Get a header consisting of `parts` joined by the '|' character.
    pub fn construct_header(&mut self, parts: [&str; N]) -> &str {
        const EXPECT_MSG: &str = "failed to write to vec for some reason";
        self.buf.clear();
        let mut writer = csv::WriterBuilder::new()
            .delimiter(DELIMITER)
            .from_writer(std::mem::take(&mut self.buf));
        writer.write_record(parts).expect(EXPECT_MSG);
        self.buf = writer.into_inner().expect(EXPECT_MSG);
        str::from_utf8(self.buf.trim_ascii_end())
            .expect("unexpected non-UTF8 made by joining UTF-8 parts")
    }
}
impl<const N: usize> CompoundHeaderReader<N> {
    pub fn new() -> Self {
        Self {
            buf: StringRecord::new(),
        }
    }
    /// Get a header consisting of the given number of parts.
    ///
    /// Fails if the number of parts is not as expected,
    /// or if the string cannot be parsed as a CSV line.
    pub fn split_header(&mut self, header: &str) -> Result<[&str; N], Error> {
        if !csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(DELIMITER)
            .from_reader(header.as_bytes())
            .read_record(&mut self.buf)?
        {
            return Err(Error::msg("empty header"));
        }
        if self.buf.len() != N {
            return Err(Error::msg(format!(
                "expected {} parts, got {}",
                N,
                self.buf.len()
            )));
        }
        Ok(std::array::from_fn::<&str, N, _>(|i| {
            self.buf.get(i).unwrap()
        }))
    }
}
