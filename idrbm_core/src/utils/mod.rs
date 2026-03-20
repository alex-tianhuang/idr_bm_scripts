mod read_to_end;
mod read_fasta;
pub use read_to_end::{read_file, read_to_end};
pub use read_fasta::read_fasta;

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

/// Workaround for box dyn stuff in a [`bumpalo::boxed::Box`].
#[macro_export]
macro_rules! box_dyn_in {
    (Box::new_in($x:expr, $arena:expr) as Box<dyn $($tokens:tt)+>) => {
        unsafe { bumpalo::boxed::Box::from_raw($arena.alloc($x) as &mut dyn $($tokens)+ as *mut dyn $($tokens)+) }
    };
}