//! Module of datatypes and associated serialization/deserialization types.
//!
//! [`sequences`] for sequence-related datatypes.
//! [`regions`] for [`RegionCsvRecord`].
mod sequences;
mod regions;
pub use sequences::{AAIndex, AAMap, AMINOACIDS, Aminoacid, FastaEntry, aa_canonical_str};
pub use regions::RegionCsvRecord;
