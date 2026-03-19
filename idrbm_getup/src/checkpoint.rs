use anyhow::{Context, Error};
use bumpalo::Bump;
use hashbrown::{DefaultHashBuilder, HashSet};
use std::{mem::ManuallyDrop, path::Path};

/// A set of headers that are already in the
/// output FASTA file.
pub struct Checkpoint<'a> {
    inner: ManuallyDrop<HashSet<&'a str, DefaultHashBuilder, &'a Bump>>,
}
impl<'a> Checkpoint<'a> {
    /// Load a checkpoint by reading the FASTA file at the given path.
    /// 
    /// If the file does not exist or is empty, return an empty checkpoint.
    pub fn load_or_empty(path: &Path, arena: &'a Bump) -> Result<Self, Error> {
        let bytes = match std::fs::read(path) {
            Ok(v) => v,
            Err(e) if matches!(e.kind(), std::io::ErrorKind::NotFound) => {
                return Ok(Self::empty(arena));
            }
            Err(e) => return Err(e.into()),
        };
        if bytes.is_empty() {
            return Ok(Self::empty(arena));
        }
        if !bytes.starts_with(b">") {
            return Err(Error::msg(format!(
                "expected fasta file starting with `>`, got `{}` character on first line",
                bytes.get(0).unwrap()
            )));
        }
        let mut headers_vec = Vec::new();
        for notification in parse_fasta_headers(&bytes) {
            let header = notification.with_context(|| {
                format!(
                    "expected UTF-8 + non-empty sequences in output FASTA file @ {}",
                    path.display()
                )
            })?;
            headers_vec.push(&*arena.alloc_str(header))
        }
        let mut headers = HashSet::new_in(arena);
        headers.extend(headers_vec);
        Ok(Self {
            inner: ManuallyDrop::new(headers),
        })
    }
    /// Given a unique pointer to a list of input IDs,
    /// remove all IDs from there that are already in this checkpoint.
    /// 
    /// The returned list re-uses the memory of the input list.
    pub fn retain<'b>(&self, ids: &'b mut [&'b str]) -> &'b [&'b str] {
        let mut start = ids.len();
        for (i, id) in ids.iter().enumerate() {
            if self.inner.contains(id) {
                start = i;
                break;
            }
        }
        if start == ids.len() {
            return ids
        }
        let mut len = start;
        let mut chunk_start = start + 1;
        for i in chunk_start..ids.len() {
            if !self.inner.contains(ids[i]) {
                continue;
            };
            if chunk_start < i {
                let chunk_range = chunk_start..i;
                ids.copy_within(chunk_range.clone(), len);
                len += chunk_range.len();
            }
            chunk_start = i + 1;
        }
        if chunk_start < ids.len() {
            let chunk_range = chunk_start..ids.len();
            ids.copy_within(chunk_range.clone(), len);
            len += chunk_range.len();
        }
        &ids[..len]
    }
    /// Make an empty checkpoint with no existing headers.
    fn empty(arena: &'a Bump) -> Self {
        Self {
            inner: ManuallyDrop::new(HashSet::new_in(arena)),
        }
    }
}
/// Iterate over a buffer with FASTA-formatted sequences
/// and only yield the headers.
///
/// Ripped from [`idrbm_core::utils::read_fasta`] module.
///
/// This function will treat the first line of this slice as a header
/// (assumes it starts with `>` character). On debug builds this crashes
/// the program if the assumption is wrong.
fn parse_fasta_headers(mut slice: &[u8]) -> impl Iterator<Item = Result<&str, Error>> {
    debug_assert!(slice.starts_with(b">"));
    std::iter::from_fn(move || {
        if slice.is_empty() {
            return None;
        }
        let working_slice: &[u8] = std::mem::replace(&mut slice, &[]);
        let Some(idx) = end_of_current_header(working_slice) else {
            return Some(Err(format!(
                "could not parse entry (empty sequence): {}",
                String::from_utf8_lossy(working_slice.strip_prefix(b">").unwrap_or(working_slice))
            )));
        };
        let (header_slice, rest) = working_slice.split_at(idx);
        debug_assert!(header_slice.starts_with(b">"));
        let header_slice = &header_slice[1..];
        let working_slice = &rest[1..];
        match next_start_of_header(working_slice) {
            Some(idx) => {
                let rest: &[u8];
                ((_, rest)) = working_slice.split_at(idx);
                slice = &rest[1..];
            }
            None => {
                slice = &[];
            }
        }
        let r = str::from_utf8(header_slice).map_err(|e| {
            format!(
                "could not parse entry ({}): {}",
                e,
                String::from_utf8_lossy(header_slice)
            )
        });
        Some(r)
    })
    .map(|r| r.map_err(Error::msg))
}
/// Find the index of the next newline from the beginning of this slice.
fn end_of_current_header(bytes: &[u8]) -> Option<usize> {
    let b = bytes.iter().find(|b| **b == b'\n')?;
    let offset = (b as *const u8 as usize) - (bytes.as_ptr() as usize);
    Some(offset)
}
/// Find the index of the next newline followed by a `>` character
/// from the beginning of this slice.
fn next_start_of_header(bytes: &[u8]) -> Option<usize> {
    let ptr = bytes.windows(2).find(|pair| *pair == b"\n>")?;
    let offset = (ptr.as_ptr() as usize) - (bytes.as_ptr() as usize);
    Some(offset)
}
