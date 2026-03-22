use bumpalo::Bump;
use hashbrown::DefaultHashBuilder;
use crate::datatypes::RegionBounds;

/// A collection of regions, grouped by protein ID.
pub type Regions<'a> = hashbrown::HashMap<&'a str, &'a [(&'a str, RegionBounds)], DefaultHashBuilder, &'a Bump>;