//! Defines [`Grouped`], [`Regions`], and [`Variants`],
//! for grouping data by protein ID without many independent
//! string allocations.
use crate::datatypes::{RegionBounds, aa_canonical_str};
use bumpalo::Bump;
use hashbrown::DefaultHashBuilder;

/// Collections of stuff, grouped by protein ID.
/// 
/// `order` must contain all keys of `mapping`.
/// Failure to do this results in runtime panics
/// when iterating over [`Grouped`].
pub struct Grouped<'a, T> {
    pub(crate) mapping: 
        hashbrown::HashMap<&'a str, T, DefaultHashBuilder, &'a Bump>,
    pub(crate) order: &'a [&'a str]
}
/// A collection of regions, grouped by protein ID.
pub type Regions<'a> = Grouped<'a, &'a [(&'a str, RegionBounds)]>;

/// A collection of variant sequences, grouped by protein ID and region ID.
pub type Variants<'a>  = Grouped<'a, &'a [(&'a str, &'a [(&'a str, &'a aa_canonical_str)])]>;

impl<'a, T: Copy> Grouped<'a, T> {
    /// Make a new [`Grouped`] from a map and order.
    /// 
    /// `order` must contain all keys of `mapping`.
    /// Failure to do this results in runtime panics
    /// when iterating over [`Grouped`].
    pub fn new(order: &'a [&'a str], mapping: hashbrown::HashMap<&'a str, T, DefaultHashBuilder, &'a Bump>) -> Self {
        assert_eq!(order.len(), mapping.len());
        Self { mapping, order }
    }
    /// Get an entry from the inner mapping.
    pub fn get(&self, protein_id: &str) -> Option<T> {
        self.mapping.get(protein_id).copied()
    }
    /// Call `contains_key` on the inner mapping.
    pub fn contains_key(&self, protein_id: &str) -> bool {
        self.mapping.contains_key(protein_id)
    }
}
impl<'a, T: Copy> IntoIterator for Grouped<'a, T> {
    type IntoIter = GroupedIter<'a, T>;
    type Item = (&'a str, T);
    fn into_iter(self) -> Self::IntoIter {
        GroupedIter {
            iter: self.order.iter(),
            mapping: self.mapping
        }
    }
}
pub struct GroupedIter<'a, T> {
    iter: std::slice::Iter<'a, &'a str>,
    mapping:
        hashbrown::HashMap<&'a str, T, DefaultHashBuilder, &'a Bump>,
}
impl<'a, T: Copy> Iterator for GroupedIter<'a, T> {
    type Item = (&'a str, T);
    fn next(&mut self) -> Option<Self::Item> {
        let protein_id = *self.iter.next()?;
        Some((protein_id, *self.mapping.get(protein_id).unwrap()))
    }
}
