//! Submodule of IDR benchmarking datatypes.
mod grouped;
mod records;
pub use grouped::{Regions, Variants, Grouped};
pub use records::{RegionBounds, RegionCsvRecord, VariantCsvRecord};