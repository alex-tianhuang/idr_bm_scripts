use bumpalo::Bump;
use hashbrown::DefaultHashBuilder;
use crate::datatypes::{RegionBounds, aa_canonical_str};

/// A collection of regions, grouped by protein ID.
pub type Regions<'a> = hashbrown::HashMap<&'a str, &'a [(&'a str, RegionBounds)], DefaultHashBuilder, &'a Bump>;
/// A collection of variant sequences, grouped by protein ID and region ID.
pub type Variants<'a> = hashbrown::HashMap<&'a str, &'a [(&'a str, &'a [(&'a str, &'a aa_canonical_str)])], DefaultHashBuilder, &'a Bump>;