//! Regression test: verify that `strata_count=None` (the default) does not
//! cause the mapper to crash when exploring SA entries with suffixes shorter
//! than the seed.
//!
//! The bug was that the FFI shim set YARA's strataCount sentinel to
//! `uint64_t(-1)` instead of `uint64_t(-1u)`, causing `getReadStrata()` to use
//! a massive strata count instead of `readSeqLength * strataRate`.  This made
//! the mapper search far more aggressively, hitting SA entries that triggered
//! an assertion (or UB with NDEBUG).

use std::io::Write;

use yara_mapper::{
    IndexerOptions, MapperOptions, ReadBatch, ReadEnd, Sensitivity, YaraIndexer, YaraMapper,
};

/// Build a FASTA with many short, highly similar contigs — the kind of
/// reference that produces SA entries with suffixes shorter than seeds.
fn write_test_fasta(path: &std::path::Path) {
    let mut f = std::fs::File::create(path).unwrap();
    // Base contig sequence (150 bp)
    let base_seq = "ATTTTATAAAGAGAAGCCTGGGGCAAAAATAAATTCAGTAATTTGTTGACTCTTCTAAAGCACATTAGTGGTGGAACTGCAACTCACCATTATTTCCTTCTAAGACCTTTGCTCTTCTCCCCAGGACTTAAGGCTCTTCAGCGTGTCTAA";
    // Create many similar contigs (vary a few positions) to produce a complex
    // FM index with boundary SA entries.
    let substitutions = [('A', 'C'), ('C', 'G'), ('G', 'T'), ('T', 'A')];
    for i in 0..200 {
        writeln!(f, ">contig_{i}").unwrap();
        let mut seq: Vec<u8> = base_seq.bytes().collect();
        // Mutate a few positions based on the contig index to create variation.
        for j in 0..3 {
            let pos = (i * 3 + j * 7) % seq.len();
            let orig = seq[pos] as char;
            for &(from, to) in &substitutions {
                if orig == from {
                    seq[pos] = to as u8;
                    break;
                }
            }
        }
        f.write_all(&seq).unwrap();
        writeln!(f).unwrap();
    }
}

#[test]
fn strata_count_sentinel_does_not_crash() {
    let tmp = tempfile::tempdir().unwrap();
    let fasta = tmp.path().join("ref.fasta");
    let prefix = tmp.path().join("ref.fasta.yara");

    // Build index.
    write_test_fasta(&fasta);
    let _indexer = YaraIndexer::build(&fasta, &prefix, &IndexerOptions::default())
        .expect("index build failed");

    // Map with Full sensitivity and default strata_count (None → sentinel).
    let options = MapperOptions {
        sensitivity: Sensitivity::Full,
        threads: 1,
        strata_count: None,
        ..Default::default()
    };
    let mapper = YaraMapper::open(&prefix, &options).expect("open failed");

    // Use the base sequence as a read pair — should map without crashing.
    let seq = b"ATTTTATAAAGAGAAGCCTGGGGCAAAAATAAATTCAGTAATTTGTTGACTCTTCTAAAGCACATTAGTGGTGGAACTGCAACTCACCATTATTTCCTTCTAAGACCTTTGCTCTTCTCCCCAGGACTTAAGGCTCTTCAGCGTGTCTAA";
    let qual = vec![b'I'; seq.len()];

    let mut batch = ReadBatch::new();
    batch.push("test_read", ReadEnd { seq, qual: &qual }, ReadEnd { seq, qual: &qual }).unwrap();

    // The key assertion: map_paired succeeds without panicking or segfaulting.
    let records = mapper.map_paired(&batch).expect("mapping failed");
    assert!(!records.is_empty(), "expected at least one alignment record");
}
