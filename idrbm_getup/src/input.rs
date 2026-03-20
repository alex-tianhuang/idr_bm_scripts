//! Module defining [`read_ids`].
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use idrbm_core::utils::{leak_vec, read_file};
use std::path::Path;
/// Read the input file of Uniprot IDs (one per line).
pub fn read_ids<'a>(path: &Path, arena: &'a Bump) -> Result<&'a mut [&'a str], Error> {
    let contents = leak_vec(read_file(path, arena)?);
    let mut lines = Vec::new_in(arena);
    for line in contents
        .split(|b| *b == b'\n')
        .filter(|slice| !slice.trim_ascii().is_empty())
    {
        let line = str::from_utf8(line.trim_ascii()).map_err(|e| {
            Error::msg(format!(
                "non-UTF8 data in input file ({}): {}",
                e,
                path.display()
            ))
        })?;
        lines.push(line)
    }
    Ok(leak_vec(lines))
}
