mod read_file;
mod read_fasta;
pub use read_file::read_file;
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