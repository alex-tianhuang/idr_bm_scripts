use crate::datatypes::{RegionBounds, aa_canonical_str};
use bumpalo::Bump;
use hashbrown::DefaultHashBuilder;

/// A collection of regions, grouped by protein ID.
pub struct Regions<'a> {
    pub(crate) mapping:
        hashbrown::HashMap<&'a str, &'a [(&'a str, RegionBounds)], DefaultHashBuilder, &'a Bump>,
    pub(crate) order: &'a [&'a str],
}
/// A collection of variant sequences, grouped by protein ID and region ID.
pub struct Variants<'a> {
    pub(crate) mapping: hashbrown::HashMap<
        &'a str,
        &'a [(&'a str, &'a [(&'a str, &'a aa_canonical_str)])],
        DefaultHashBuilder,
        &'a Bump,
    >,
    pub(crate) order: &'a [&'a str],
}
impl<'a> Regions<'a> {
    /// Get an entry from the inner mapping.
    pub fn get(&self, protein_id: &str) -> Option<&'a [(&'a str, RegionBounds)]> {
        self.mapping.get(protein_id).copied()
    }
}
impl<'a> IntoIterator for Regions<'a> {
    type IntoIter = OrderedMapIter<'a, &'a [(&'a str, RegionBounds)]>;
    type Item = (&'a str, &'a [(&'a str, RegionBounds)]);
    fn into_iter(self) -> Self::IntoIter {
        OrderedMapIter {
            iter: self.order.iter(),
            mapping: self.mapping
        }
    }
}
impl<'a> IntoIterator for Variants<'a> {
    type IntoIter = OrderedMapIter<'a, &'a [(&'a str, &'a [(&'a str, &'a aa_canonical_str)])]>;
    type Item = (&'a str, &'a [(&'a str, &'a [(&'a str, &'a aa_canonical_str)])]);
    fn into_iter(self) -> Self::IntoIter {
        OrderedMapIter {
            iter: self.order.iter(),
            mapping: self.mapping
        }
    }
}
pub struct OrderedMapIter<'a, T> {
    iter: std::slice::Iter<'a, &'a str>,
    mapping:
        hashbrown::HashMap<&'a str, T, DefaultHashBuilder, &'a Bump>,
}
impl<'a, T: Copy> Iterator for OrderedMapIter<'a, T> {
    type Item = (&'a str, T);
    fn next(&mut self) -> Option<Self::Item> {
        let protein_id = *self.iter.next()?;
        Some((protein_id, *self.mapping.get(protein_id).unwrap()))
    }
}
