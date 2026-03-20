//! Two submodules for communicating with Uniprot web APIs.
mod submit;
mod results;
pub use submit::submit;
pub use results::{results, alloc_error_handler};