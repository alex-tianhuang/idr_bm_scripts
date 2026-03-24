//! Some utilities with extremely obtuse names
//! for the csv module.
use crate::{csv::duplicate_rule::DuplicateRule, datatypes::{Grouped, RegionMap, VariantMap}};
use bumpalo::{Bump, collections::Vec};
use hashbrown::{DefaultHashBuilder, HashMap};
use std::fmt::Debug;

pub fn fold_l1<'a, 'b, V, T: Debug>(
    protein_id: &str,
    item: V,
    default_items: impl Fn() -> T,
    fold_items: impl Fn(&mut T, V),
    local_mapping: &mut HashMap<&'a str, T, DefaultHashBuilder, &'b Bump>,
    local_order: &mut Vec<'b, &'a str>,
    arena: &'a Bump,
) {
    let protein_group = match local_mapping.get_mut(protein_id) {
        Some(group) => group,
        None => {
            let protein_id = &*arena.alloc_str(protein_id);
            local_order.push(protein_id);
            local_mapping.try_insert(protein_id, default_items()).unwrap()
        }
    };
    fold_items(protein_group, item)
}
pub fn fold_l2<'a, 'b, V, T: Debug>(
    protein_id: &str,
    region_id: &str,
    item: V,
    default_items: impl Fn() -> T,
    fold_items: impl Fn(&mut T, V),
    local_mapping: &mut HashMap<&'a str, Vec<'b, (&'a str, T)>, DefaultHashBuilder, &'b Bump>,
    local_order: &mut Vec<'b, &'a str>,
    local_arena: &'b Bump,
    arena: &'a Bump,
) {
    fold_l1(
        protein_id,
        (region_id, item),
        || Vec::new_in(local_arena),
        |protein_group, (region_id, item)| {
            let (_, region_group) =
                match protein_group.iter_mut().find(|(rid, _)| *rid == region_id) {
                    Some(tuple) => tuple,
                    None => {
                        let region_id = &*arena.alloc_str(region_id);
                        protein_group.push((region_id, default_items()));
                        protein_group.last_mut().unwrap()
                    }
                };
            fold_items(region_group, item)
        },
        local_mapping,
        local_order,
        arena,
    );
}
pub fn collect_l2<'a, 'b, T: Debug>(
    protein_id: &str,
    region_id: &str,
    item: T,
    duplicate_rule: DuplicateRule,
    local_mapping: &mut HashMap<&'a str, Vec<'b, (&'a str, T)>, DefaultHashBuilder, &'b Bump>,
    local_order: &mut Vec<'b, &'a str>,
    local_arena: &'b Bump,
    arena: &'a Bump,
) {
    fold_l1(
        protein_id,
        (region_id, item),
        || Vec::new_in(local_arena),
        |protein_group, (region_id, item)| {
            match protein_group.iter_mut().find(|(rid, _)| *rid == region_id) {
                Some((_, slot)) => {
                    if matches!(duplicate_rule, DuplicateRule::LastWins) {
                        *slot = item;
                    }
                },
                None => {
                    let region_id = &*arena.alloc_str(region_id);
                    protein_group.push((region_id, item));
                }
            };
        },
        local_mapping,
        local_order,
        arena,
    );
}
pub fn collect_l3<'a, 'b, T: Debug>(
    protein_id: &str,
    region_id: &str,
    variant_id: &str,
    item: T,
    duplicate_rule: DuplicateRule,
    local_mapping: &mut HashMap<&'a str, Vec<'b, (&'a str, Vec<'b, (&'a str, T)>)>, DefaultHashBuilder, &'b Bump>,
    local_order: &mut Vec<'b, &'a str>,
    local_arena: &'b Bump,
    arena: &'a Bump,
) {
    fold_l2(
        protein_id,
        region_id,
        (variant_id, item),
        || Vec::new_in(local_arena),
        |region_group, (variant_id, item)| {
            match region_group.iter_mut().find(|(vid, _)| *vid == variant_id) {
                Some((_, slot)) => {
                    if matches!(duplicate_rule, DuplicateRule::LastWins) {
                        *slot = item;
                    }
                },
                None => {
                    let variant_id = &*arena.alloc_str(variant_id);
                    region_group.push((variant_id, item));
                }
            };
        },
        local_mapping,
        local_order,
        local_arena,
        arena,
    );
}
pub fn finish_l2<'a, T, R>(
    finish_t: impl Fn(&T) -> R,
    local_mapping: HashMap<&'a str, Vec<'_, (&'a str, T)>, DefaultHashBuilder, &Bump>,
    local_order: Vec<'_, &'a str>,
    arena: &'a Bump,
) -> RegionMap<'a, R> {
    let mut mapping = HashMap::with_capacity_in(local_mapping.len(), arena);
    for (protein_id, protein_group) in local_mapping.iter() {
        // SAFETY: keys are unique because they come from a map
        unsafe {
            mapping.insert_unique_unchecked(
                *protein_id,
                &*arena.alloc_slice_fill_iter(
                    protein_group
                        .iter()
                        .map(|(region_id, region_group)| (*region_id, finish_t(region_group))),
                ),
            )
        };
    }
    let ret = Grouped {
        mapping,
        order: arena.alloc_slice_copy(&local_order),
    };
    std::mem::forget(local_mapping);
    std::mem::forget(local_order);
    ret
}
pub fn finish_l3<'a, T, R>(
    finish_t: impl Fn(&T) -> R,
    local_mapping: HashMap<
        &'a str,
        Vec<'_, (&'a str, Vec<'_, (&'a str, T)>)>,
        DefaultHashBuilder,
        &Bump,
    >,
    local_order: Vec<'_, &'a str>,
    arena: &'a Bump,
) -> VariantMap<'a, R> {
    finish_l2(
        |region_group| {
            &*arena.alloc_slice_fill_iter(
                region_group
                    .iter()
                    .map(|(variant_id, variant_group)| (*variant_id, finish_t(variant_group))),
            )
        },
        local_mapping,
        local_order,
        arena,
    )
}
