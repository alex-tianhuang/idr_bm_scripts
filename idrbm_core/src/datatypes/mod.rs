//! Module of datatypes and associated serialization/deserialization types.
//!
//! [`sequences`] for sequence-related datatypes.
//! [`idrbm`] for IDR benchmarking datatypes specifically.
//! [`grouped`] for collections of records.
mod sequences;
mod idrbm;
pub use sequences::{AAIndex, AAMap, AMINOACIDS, Aminoacid, FastaEntry, aa_canonical_str};
pub use idrbm::{RegionBounds, RegionMap, VariantMap, Grouped};