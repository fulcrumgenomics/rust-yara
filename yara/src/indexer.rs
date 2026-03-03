use std::ffi::CStr;
use std::path::{Path, PathBuf};
use std::ptr::NonNull;

use crate::error::YaraError;
use crate::ffi_helpers::{collect_contig_lengths, collect_contig_names, path_to_cstring};

/// Options for building a YARA FM index.
#[derive(Debug, Clone, Default)]
pub struct IndexerOptions {
    /// Temporary directory for intermediate files during index construction.
    /// `None` uses the output directory's parent.
    pub tmp_dir: Option<PathBuf>,
    /// Verbosity level (0 = silent, 1 = progress, 2 = debug).
    pub verbose: u32,
}

/// Handle to a built YARA FM index.
///
/// After calling [`YaraIndexer::build`], the on-disk index files have been
/// written and the handle retains contig metadata for querying.
///
/// # Thread safety
///
/// [`YaraIndexer`] is [`Send`] but not [`Sync`], consistent with
/// [`crate::YaraMapper`].  Note that [`YaraIndexer::build`] modifies
/// the `TMPDIR` environment variable (restored via RAII), so concurrent
/// `build` calls from multiple threads are not safe.
pub struct YaraIndexer {
    handle: NonNull<yara_mapper_sys::YaraIndexerHandle>,
}

// SAFETY: The C++ handle owns all its memory and can be moved between threads.
// Post-build accessors are read-only, but we keep !Sync for consistency with
// YaraMapper and because build() modifies global state (TMPDIR).
unsafe impl Send for YaraIndexer {}

impl YaraIndexer {
    /// Build a YARA FM index from a FASTA file.
    ///
    /// `fasta_path` is the path to the input FASTA/FASTQ reference file.
    /// `output_prefix` is the file prefix for index output (e.g., `ref/genome`
    /// produces `ref/genome.txt`, `ref/genome.rid`, `ref/genome.txt.size`,
    /// `ref/genome.sa`, `ref/genome.lf`).
    ///
    /// # Errors
    ///
    /// Returns [`YaraError::IndexBuild`] if the FASTA file cannot be opened,
    /// if the paths contain non-UTF-8 characters or null bytes, or if the
    /// underlying C++ indexer fails (e.g., insufficient memory or disk space).
    pub fn build<P: AsRef<Path>>(
        fasta_path: P,
        output_prefix: P,
        options: &IndexerOptions,
    ) -> Result<Self, YaraError> {
        let fasta_cstr = path_to_cstring(fasta_path.as_ref(), "fasta_path", YaraError::IndexBuild)?;
        let prefix_cstr =
            path_to_cstring(output_prefix.as_ref(), "output_prefix", YaraError::IndexBuild)?;
        let tmp_dir_cstr = options
            .tmp_dir
            .as_ref()
            .map(|p| path_to_cstring(p, "tmp_dir", YaraError::IndexBuild))
            .transpose()?;

        let ffi_opts = yara_mapper_sys::YaraIndexerOptions {
            output_prefix: prefix_cstr.as_ptr(),
            tmp_dir: tmp_dir_cstr.as_ref().map_or(std::ptr::null(), |c| c.as_ptr()),
            #[expect(clippy::cast_possible_wrap, reason = "FFI boundary: verbose fits in i32")]
            verbose: options.verbose as std::os::raw::c_int,
        };

        let mut error_buf = vec![0u8; 1024];

        let handle = unsafe {
            yara_mapper_sys::yara_indexer_build(
                fasta_cstr.as_ptr(),
                &ffi_opts,
                error_buf.as_mut_ptr().cast(),
                error_buf.len(),
            )
        };

        NonNull::new(handle).map(|h| Self { handle: h }).ok_or_else(|| {
            let msg = unsafe { CStr::from_ptr(error_buf.as_ptr().cast()) };
            YaraError::IndexBuild(msg.to_string_lossy().into_owned())
        })
    }

    /// Number of contigs in the built index.
    #[must_use]
    pub fn contig_count(&self) -> usize {
        unsafe { yara_mapper_sys::yara_indexer_contig_count(self.handle.as_ptr()) }
    }

    /// Contig names from the indexed reference.
    #[must_use]
    pub fn contig_names(&self) -> Vec<String> {
        let n = self.contig_count();
        unsafe {
            collect_contig_names(n, |i| {
                yara_mapper_sys::yara_indexer_contig_name(self.handle.as_ptr(), i)
            })
        }
    }

    /// Contig lengths from the indexed reference.
    #[must_use]
    pub fn contig_lengths(&self) -> Vec<usize> {
        let n = self.contig_count();
        collect_contig_lengths(n, |i| unsafe {
            yara_mapper_sys::yara_indexer_contig_length(self.handle.as_ptr(), i)
        })
    }
}

impl Drop for YaraIndexer {
    fn drop(&mut self) {
        unsafe { yara_mapper_sys::yara_indexer_close(self.handle.as_ptr()) }
    }
}
