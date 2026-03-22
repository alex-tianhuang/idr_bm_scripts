//! Module of datatypes and associated serialization/deserialization types.
//!
//! [`sequences`] for sequence-related datatypes.
//! [`records`] for CSV record types for one sequence or region.
//! [`grouped`] for collections of records.
mod sequences;
mod records;
mod grouped;
pub use sequences::{AAIndex, AAMap, AMINOACIDS, Aminoacid, FastaEntry, aa_canonical_str};
pub use records::{RegionCsvRecord, VariantCsvRecord, RegionBounds};
pub use grouped::{Regions, Variants};