use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::path::Path;

use crate::error::YaraError;

/// Convert a [`Path`] to a [`CString`] for passing across the FFI boundary.
///
/// # Errors
///
/// Returns an error via `make_err` if the path contains non-UTF-8 characters
/// or embedded null bytes.
pub(crate) fn path_to_cstring(
    path: &Path,
    label: &str,
    make_err: fn(String) -> YaraError,
) -> Result<CString, YaraError> {
    let s = path.to_str().ok_or_else(|| make_err(format!("{label} is not valid UTF-8")))?;
    CString::new(s).map_err(|e| make_err(format!("{label} contains null byte: {e}")))
}

/// Collect contig names from an FFI handle by calling `name_fn` for each index.
///
/// # Safety
///
/// `name_fn` must return either a null pointer or a valid, null-terminated C
/// string whose lifetime extends at least until this function returns.
pub(crate) unsafe fn collect_contig_names(
    count: usize,
    name_fn: impl Fn(usize) -> *const c_char,
) -> Vec<String> {
    (0..count)
        .map(|i| {
            let ptr = name_fn(i);
            if ptr.is_null() {
                String::new()
            } else {
                unsafe { CStr::from_ptr(ptr) }.to_string_lossy().into_owned()
            }
        })
        .collect()
}

/// Convert a byte slice of DNA bases or quality scores to a [`CString`].
///
/// # Safety
///
/// The caller must ensure the bytes do not contain interior null bytes.
/// This is guaranteed for ASCII DNA bases (ACGTN, a-t) and phred+33
/// quality values (0x21..0x7E).
pub(crate) fn bytes_to_cstring(bytes: &[u8]) -> CString {
    unsafe { CString::from_vec_unchecked(bytes.to_vec()) }
}

/// Collect contig lengths from an FFI handle by calling `length_fn` for each index.
pub(crate) fn collect_contig_lengths(
    count: usize,
    length_fn: impl Fn(usize) -> usize,
) -> Vec<usize> {
    (0..count).map(length_fn).collect()
}
