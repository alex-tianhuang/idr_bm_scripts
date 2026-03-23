//! Submodule of IDR benchmarking datatypes.
mod grouped;
mod regions;
pub use regions::RegionBounds;
pub use grouped::{Grouped, RegionMap, VariantMap};
