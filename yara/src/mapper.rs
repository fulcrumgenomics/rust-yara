use std::ffi::{CStr, CString};
use std::path::Path;
use std::ptr::NonNull;

use crate::error::YaraError;
use crate::ffi_helpers::{
    bytes_to_cstring, collect_contig_lengths, collect_contig_names, path_to_cstring,
};
use crate::options::{MapperOptions, SecondaryMode};
use crate::record::{CigarOp, YaraRecord};

/// One end (read 1 or read 2) of a paired-end read.
///
/// Sequence should be ASCII DNA (ACGTN) and quality should be phred+33.
#[derive(Debug, Clone, Copy)]
pub struct ReadEnd<'a> {
    /// DNA sequence bytes (ASCII).
    pub seq: &'a [u8],
    /// Base quality bytes (phred+33 ASCII).
    pub qual: &'a [u8],
}

/// A batch of paired-end reads to map.
///
/// All slices must have the same length (`count`).  Sequences are ASCII
/// DNA strings (ACGTN) and qualities are phred+33 ASCII strings.
#[derive(Default)]
pub struct ReadBatch {
    names: Vec<CString>,
    r1_seqs: Vec<CString>,
    r1_quals: Vec<CString>,
    r2_seqs: Vec<CString>,
    r2_quals: Vec<CString>,
}

impl ReadBatch {
    /// Create a new empty batch.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a batch with pre-allocated capacity for `n` pairs.
    #[must_use]
    pub fn with_capacity(n: usize) -> Self {
        Self {
            names: Vec::with_capacity(n),
            r1_seqs: Vec::with_capacity(n),
            r1_quals: Vec::with_capacity(n),
            r2_seqs: Vec::with_capacity(n),
            r2_quals: Vec::with_capacity(n),
        }
    }

    /// Add a read pair to the batch.
    ///
    /// # Errors
    ///
    /// Returns [`YaraError::InvalidInput`] if any input contains a null byte
    /// or if sequence and quality lengths do not match.
    pub fn push(&mut self, name: &str, r1: ReadEnd<'_>, r2: ReadEnd<'_>) -> Result<(), YaraError> {
        if r1.seq.len() != r1.qual.len() {
            return Err(YaraError::InvalidInput(format!(
                "r1 seq/qual length mismatch: {} vs {}",
                r1.seq.len(),
                r1.qual.len()
            )));
        }
        if r2.seq.len() != r2.qual.len() {
            return Err(YaraError::InvalidInput(format!(
                "r2 seq/qual length mismatch: {} vs {}",
                r2.seq.len(),
                r2.qual.len()
            )));
        }
        self.names
            .push(CString::new(name).map_err(|e| YaraError::InvalidInput(format!("name: {e}")))?);
        self.r1_seqs.push(bytes_to_cstring(r1.seq));
        self.r1_quals.push(bytes_to_cstring(r1.qual));
        self.r2_seqs.push(bytes_to_cstring(r2.seq));
        self.r2_quals.push(bytes_to_cstring(r2.qual));
        Ok(())
    }

    /// Number of read pairs in the batch.
    #[must_use]
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Whether the batch is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// Clear all reads, retaining allocated capacity for reuse.
    pub fn clear(&mut self) {
        self.names.clear();
        self.r1_seqs.clear();
        self.r1_quals.clear();
        self.r2_seqs.clear();
        self.r2_quals.clear();
    }
}

/// Handle to a loaded YARA mapper with a pre-built FM index.
///
/// The mapper is configured at construction time and can be used to map
/// multiple batches of reads.  The FM index stays loaded in memory for
/// the lifetime of this handle.
///
/// # Thread safety
///
/// The underlying C++ mapper uses OpenMP internally but is not safe for
/// concurrent `map_paired` calls.  [`YaraMapper`] is [`Send`] but not
/// [`Sync`].
pub struct YaraMapper {
    handle: NonNull<yara_mapper_sys::YaraMapperHandle>,
    secondary_mode: SecondaryMode,
}

// SAFETY: The C++ handle owns all its memory and can be moved between threads.
// It is NOT safe to call map_paired concurrently (OpenMP parallelism is
// internal to each call), so we implement Send but not Sync.
unsafe impl Send for YaraMapper {}

impl YaraMapper {
    /// Open a pre-built YARA index and create a mapper.
    ///
    /// `index_prefix` is the path prefix used when building the index (e.g.,
    /// `ref/hla.fasta` if the index files are `ref/hla.fasta.yara.*`).
    ///
    /// # Errors
    ///
    /// Returns [`YaraError::IndexOpen`] if the index files cannot be found or
    /// loaded, or if `index_prefix` contains non-UTF-8 characters or null bytes.
    pub fn open<P: AsRef<Path>>(
        index_prefix: P,
        options: &MapperOptions,
    ) -> Result<Self, YaraError> {
        let prefix_cstr =
            path_to_cstring(index_prefix.as_ref(), "index_prefix", YaraError::IndexOpen)?;

        let ffi_opts = options.to_ffi();
        let mut error_buf = vec![0u8; 1024];

        let handle = unsafe {
            yara_mapper_sys::yara_mapper_open(
                prefix_cstr.as_ptr(),
                &ffi_opts,
                error_buf.as_mut_ptr().cast(),
                error_buf.len(),
            )
        };

        NonNull::new(handle)
            .map(|h| Self { handle: h, secondary_mode: options.secondary_mode })
            .ok_or_else(|| {
                let msg = unsafe { CStr::from_ptr(error_buf.as_ptr().cast()) };
                YaraError::IndexOpen(msg.to_string_lossy().into_owned())
            })
    }

    /// Map a batch of paired-end reads and return alignment records.
    ///
    /// The returned records include primary and (depending on options)
    /// secondary alignments, plus unmapped entries for reads that could
    /// not be aligned.
    ///
    /// # Errors
    ///
    /// Returns [`YaraError::Mapping`] if the underlying C++ mapper encounters
    /// an error during alignment.
    pub fn map_paired(&self, reads: &ReadBatch) -> Result<Vec<YaraRecord>, YaraError> {
        if reads.is_empty() {
            return Ok(Vec::new());
        }

        // Build pointer arrays for the C API.
        let name_ptrs = cstring_ptrs(&reads.names);
        let r1_seq_ptrs = cstring_ptrs(&reads.r1_seqs);
        let r1_qual_ptrs = cstring_ptrs(&reads.r1_quals);
        let r2_seq_ptrs = cstring_ptrs(&reads.r2_seqs);
        let r2_qual_ptrs = cstring_ptrs(&reads.r2_quals);

        let batch = yara_mapper_sys::YaraReadBatch {
            names: name_ptrs.as_ptr(),
            r1_seqs: r1_seq_ptrs.as_ptr(),
            r1_quals: r1_qual_ptrs.as_ptr(),
            r2_seqs: r2_seq_ptrs.as_ptr(),
            r2_quals: r2_qual_ptrs.as_ptr(),
            count: reads.len(),
        };

        let capacity = reads.len() * records_per_pair(self.secondary_mode);
        // SAFETY: YaraAlignmentRecord is #[repr(C)] POD — all-zeros is a valid
        // representation (null pointers, zero scalars).  The C++ shim also
        // memsets the buffer, but we zero here defensively.
        let mut out_records: Vec<yara_mapper_sys::YaraAlignmentRecord> =
            vec![unsafe { std::mem::zeroed() }; capacity];

        let mut error_buf = [0u8; 1024];

        let count = unsafe {
            yara_mapper_sys::yara_mapper_map_paired(
                self.handle.as_ptr(),
                &batch,
                out_records.as_mut_ptr(),
                capacity,
                error_buf.as_mut_ptr().cast(),
                error_buf.len(),
            )
        };

        if count < 0 {
            // C++ already freed any partially written records before returning -1.
            let msg = unsafe { CStr::from_ptr(error_buf.as_ptr().cast()) };
            return Err(YaraError::Mapping(msg.to_string_lossy().into_owned()));
        }

        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "count is non-negative and bounded by capacity (which is usize)"
        )]
        let n = count as usize;

        // Convert each C record to an owned Rust type, then batch-free the
        // C++ memory.  Pool-managed fields (seq/qual/cigar) are freed when the
        // pool is cleared on the next mapPaired call; only XA strings are freed
        // by free_records.
        let results: Vec<YaraRecord> = out_records[..n].iter().map(convert_record).collect();

        unsafe {
            yara_mapper_sys::yara_mapper_free_records(out_records.as_mut_ptr(), n);
        }

        Ok(results)
    }

    /// Number of reference contigs in the loaded index.
    #[must_use]
    pub fn contig_count(&self) -> usize {
        unsafe { yara_mapper_sys::yara_mapper_contig_count(self.handle.as_ptr()) }
    }

    /// Reference contig names (for SAM/BAM header construction).
    #[must_use]
    pub fn contig_names(&self) -> Vec<String> {
        let n = self.contig_count();
        unsafe {
            collect_contig_names(n, |i| {
                yara_mapper_sys::yara_mapper_contig_name(self.handle.as_ptr(), i)
            })
        }
    }

    /// Reference contig lengths.
    #[must_use]
    pub fn contig_lengths(&self) -> Vec<usize> {
        let n = self.contig_count();
        collect_contig_lengths(n, |i| unsafe {
            yara_mapper_sys::yara_mapper_contig_length(self.handle.as_ptr(), i)
        })
    }
}

impl Drop for YaraMapper {
    fn drop(&mut self) {
        unsafe { yara_mapper_sys::yara_mapper_close(self.handle.as_ptr()) }
    }
}

/// Upper bound on records per read pair, depending on secondary mode.
/// `Tag`/`Omit`: each read end produces exactly one record, so 2 per pair.
/// `Record`: secondaries are separate records; 10 per pair is generous.
fn records_per_pair(mode: SecondaryMode) -> usize {
    match mode {
        SecondaryMode::Tag | SecondaryMode::Omit => 2,
        SecondaryMode::Record => 10,
    }
}

/// Collect raw pointers from a slice of `CString` values for the C API.
fn cstring_ptrs(strings: &[CString]) -> Vec<*const i8> {
    strings.iter().map(|s| s.as_ptr()).collect()
}

/// Convert a single FFI record to an owned Rust record.
fn convert_record(rec: &yara_mapper_sys::YaraAlignmentRecord) -> YaraRecord {
    // CIGAR
    let cigar = if !rec.cigar.is_null() && rec.cigar_len > 0 {
        let slice = unsafe { std::slice::from_raw_parts(rec.cigar, rec.cigar_len as usize) };
        slice.iter().map(|&encoded| CigarOp::from_bam(encoded)).collect()
    } else {
        Vec::new()
    };

    // Sequence
    let seq = if !rec.seq.is_null() && rec.seq_len > 0 {
        let slice =
            unsafe { std::slice::from_raw_parts(rec.seq.cast::<u8>(), rec.seq_len as usize) };
        Some(slice.to_vec())
    } else {
        None
    };

    // Quality
    let qual = if !rec.qual.is_null() && rec.seq_len > 0 {
        let slice =
            unsafe { std::slice::from_raw_parts(rec.qual.cast::<u8>(), rec.seq_len as usize) };
        Some(slice.to_vec())
    } else {
        None
    };

    // XA tag
    let xa = if rec.xa.is_null() {
        None
    } else {
        Some(unsafe { CStr::from_ptr(rec.xa) }.to_string_lossy().into_owned())
    };

    YaraRecord {
        read_pair_index: rec.read_pair_index,
        is_read1: rec.is_read1 != 0,
        contig_id: rec.contig_id,
        pos: rec.pos,
        is_reverse: rec.is_reverse != 0,
        is_secondary: rec.is_secondary != 0,
        is_unmapped: rec.is_unmapped != 0,
        mapq: rec.mapq,
        nm: rec.nm,
        x0: rec.x0,
        x1: rec.x1,
        mate_contig_id: rec.mate_contig_id,
        mate_pos: rec.mate_pos,
        tlen: rec.tlen,
        flag: rec.flag,
        cigar,
        seq,
        qual,
        xa,
    }
}
