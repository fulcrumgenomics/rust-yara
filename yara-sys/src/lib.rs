//! Low-level FFI bindings to the YARA read mapper C++ shim.
//!
//! This crate compiles a C++ shim layer that wraps YARA's SeqAn2-based template code
//! and exposes a C-compatible API for loading FM indices, ingesting reads, and retrieving
//! alignment results as flat structs.

#![allow(non_camel_case_types)]

use std::os::raw::{c_char, c_int};

// ---------------------------------------------------------------------------
// Opaque handle
// ---------------------------------------------------------------------------

/// Opaque handle to a fully-configured YARA Mapper instance.
/// Allocated and owned by the C++ side; freed via [`yara_mapper_close`].
#[repr(C)]
pub struct YaraMapperHandle {
    _opaque: [u8; 0],
}

// ---------------------------------------------------------------------------
// Option / input / output structs (must match yara_shim.h exactly)
// ---------------------------------------------------------------------------

/// Mapper configuration passed to [`yara_mapper_open`].
#[repr(C)]
#[derive(Debug, Clone)]
pub struct YaraMapperOptions {
    pub error_rate: f32,
    pub strata_rate: f32,
    pub strata_count: i32,
    pub sensitivity: i32,
    pub threads: i32,
    pub secondary_mode: i32,
    pub align_secondary: i32,
    pub verify_matches: i32,
    pub verbose: i32,
    pub library_length: i32,
    pub library_dev: i32,
}

/// A batch of paired-end reads to map.
#[repr(C)]
#[derive(Debug)]
pub struct YaraReadBatch {
    pub names: *const *const c_char,
    pub r1_seqs: *const *const c_char,
    pub r1_quals: *const *const c_char,
    pub r2_seqs: *const *const c_char,
    pub r2_quals: *const *const c_char,
    pub count: usize,
}

/// A single alignment record returned by the mapper.
#[repr(C)]
#[derive(Debug, Clone)]
#[allow(clippy::pub_underscore_fields)]
pub struct YaraAlignmentRecord {
    // Read identification
    pub read_pair_index: u32,
    pub is_read1: u8,

    // Alignment position
    pub contig_id: u32,
    pub pos: u32,
    pub is_reverse: u8,
    pub is_secondary: u8,
    pub is_unmapped: u8,

    // Alignment quality
    pub mapq: u8,
    pub nm: u8,
    pub x0: u16,
    pub x1: u16,

    // Mate info
    pub mate_contig_id: u32,
    pub mate_pos: u32,
    pub tlen: i32,
    pub flag: u16,

    // CIGAR (BAM-encoded uint32 array: op | len<<4)
    pub cigar: *mut u32,
    pub cigar_len: u32,

    // Sequence and quality (NULL for secondaries)
    pub seq: *const c_char,
    pub qual: *const c_char,
    pub seq_len: u32,

    // XA tag string (NULL if none)
    pub xa: *const c_char,

    // Internal: if set, seq/qual/cigar are pool-managed (not individually freed).
    pub _pool_managed: u8,
}

// ---------------------------------------------------------------------------
// FFI function declarations
// ---------------------------------------------------------------------------

// SAFETY: These functions are implemented in C++ (yara_shim.cpp) and compiled
// via build.rs. The caller is responsible for passing valid pointers and
// respecting the lifetime contracts documented in yara_shim.h.
unsafe extern "C" {
    /// Create a mapper and load the FM index.  Returns null on error, with
    /// `error_buf` populated.
    pub fn yara_mapper_open(
        index_prefix: *const c_char,
        opts: *const YaraMapperOptions,
        error_buf: *mut c_char,
        error_buf_len: usize,
    ) -> *mut YaraMapperHandle;

    /// Map a batch of paired-end reads.  Returns the number of output records
    /// written to `out`, or -1 on error.
    pub fn yara_mapper_map_paired(
        handle: *mut YaraMapperHandle,
        reads: *const YaraReadBatch,
        out: *mut YaraAlignmentRecord,
        out_capacity: usize,
        error_buf: *mut c_char,
        error_buf_len: usize,
    ) -> i64;

    /// Number of reference contigs in the loaded index.
    pub fn yara_mapper_contig_count(handle: *const YaraMapperHandle) -> usize;

    /// Name of the contig at `idx`. The returned pointer is valid for the
    /// lifetime of the handle.
    pub fn yara_mapper_contig_name(handle: *const YaraMapperHandle, idx: usize) -> *const c_char;

    /// Length of the contig at `idx`.
    pub fn yara_mapper_contig_length(handle: *const YaraMapperHandle, idx: usize) -> usize;

    /// Free the mapper and all associated memory.
    pub fn yara_mapper_close(handle: *mut YaraMapperHandle);

    /// Free a record batch's C++-owned memory (cigar arrays, seq/qual strings,
    /// XA strings).  Must be called after processing each batch.
    pub fn yara_mapper_free_records(records: *mut YaraAlignmentRecord, count: usize);

    /// Free a single record's C++-owned memory (cigar, seq, qual, xa).
    pub fn yara_mapper_free_record(record: *mut YaraAlignmentRecord);
}

// ---------------------------------------------------------------------------
// Indexer: opaque handle
// ---------------------------------------------------------------------------

/// Opaque handle to a YARA Indexer instance.
/// Allocated and owned by the C++ side; freed via [`yara_indexer_close`].
#[repr(C)]
pub struct YaraIndexerHandle {
    _opaque: [u8; 0],
}

// ---------------------------------------------------------------------------
// Indexer: options struct (must match yara_indexer_shim.h exactly)
// ---------------------------------------------------------------------------

/// Indexer configuration passed to [`yara_indexer_build`].
#[repr(C)]
#[derive(Debug)]
pub struct YaraIndexerOptions {
    pub output_prefix: *const c_char,
    pub tmp_dir: *const c_char,
    pub verbose: c_int,
}

// ---------------------------------------------------------------------------
// Indexer: FFI function declarations
// ---------------------------------------------------------------------------

// SAFETY: These functions are implemented in C++ (yara_indexer_shim.cpp) and
// compiled via build.rs. The caller is responsible for passing valid pointers
// and respecting the lifetime contracts documented in yara_indexer_shim.h.
unsafe extern "C" {
    /// Build an FM index from a FASTA file. Returns null on error, with
    /// `error_buf` populated.
    pub fn yara_indexer_build(
        fasta_path: *const c_char,
        opts: *const YaraIndexerOptions,
        error_buf: *mut c_char,
        error_buf_len: usize,
    ) -> *mut YaraIndexerHandle;

    /// Number of contigs in the built index.
    pub fn yara_indexer_contig_count(handle: *const YaraIndexerHandle) -> usize;

    /// Name of the contig at `idx`. The returned pointer is valid for the
    /// lifetime of the handle.
    pub fn yara_indexer_contig_name(handle: *const YaraIndexerHandle, idx: usize) -> *const c_char;

    /// Length of the contig at `idx`.
    pub fn yara_indexer_contig_length(handle: *const YaraIndexerHandle, idx: usize) -> usize;

    /// Free the indexer handle and all associated memory.
    pub fn yara_indexer_close(handle: *mut YaraIndexerHandle);
}
