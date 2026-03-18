//! Module of datatypes and associated serialization/deserialization types.
//!
//! [`sequences`] for sequence-related datatypes.
mod sequences;
pub use sequences::{AAIndex, AAMap, AMINOACIDS, Aminoacid, FastaEntry, aa_canonical_str};
