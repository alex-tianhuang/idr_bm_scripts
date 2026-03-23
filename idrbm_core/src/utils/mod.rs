mod read_fasta;
mod read_to_end;
mod read_csvs;
mod compound_header;
pub use compound_header::{CompoundHeaderWriter, CompoundHeaderReader};
pub use read_csvs::{read_regions, read_variants, read_region_csv_template, read_variant_csv_template, assert_n_cols};
pub use read_fasta::{read_fasta, read_fasta_to_map};
pub use read_to_end::{read_file, read_to_end};

/// Leak a [`bumpalo::collections::Vec`] managed by an arena into a slice
/// that lives for as long as the arena does not reset.
///
/// Analogue of [`std::vec::Vec::leak`].
pub fn leak_vec<'a, T>(buf: bumpalo::collections::Vec<'a, T>) -> &'a mut [T] {
    let mut buf = std::mem::ManuallyDrop::new(buf);
    // SAFETY: the arena must de-allocate before this vec does.
    //         Also, do to the use of `ManuallyDrop`, these bytes
    //         are not de-allocated and reused by the arena.
    unsafe { std::mem::transmute::<&mut [T], &'a mut [T]>(buf.as_mut_slice()) }
}

/// Leak a [`bumpalo::collections::String`] managed by an arena into a `str`
/// that lives for as long as the arena does not reset.
///
/// Analogue of [`std::string::String::leak`].
pub fn leak_str<'a>(s: bumpalo::collections::String<'a>) -> &'a mut str {
    let mut s = std::mem::ManuallyDrop::new(s);
    // SAFETY: the arena must de-allocate before this vec does.
    //         Also, do to the use of `ManuallyDrop`, these bytes
    //         are not de-allocated and reused by the arena.
    unsafe { std::mem::transmute::<&mut str, &'a mut str>(s.as_mut_str()) }
}

/// Workaround for box dyn stuff in a [`bumpalo::boxed::Box`].
#[macro_export]
macro_rules! box_dyn_in {
    (Box::new_in($x:expr, $arena:expr) as Box<dyn $trait_:ident>) => {
        unsafe {
            bumpalo::boxed::Box::from_raw($arena.alloc($x) as &mut dyn $trait_ as *mut dyn $trait_)
        }
    };
}
