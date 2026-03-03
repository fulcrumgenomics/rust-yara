# Yara Rust Port Implementation Plan

This document outlines a plan for porting Yara to Rust, focused specifically on HLA calling workflows (e.g., OptiType integration).

## Overview

**Goal:** Create a Rust implementation of Yara's core read alignment functionality optimized for HLA typing.

**Scope:** Single-end and paired-end mapping to small reference genomes (HLA ~6-7MB)

**Out of Scope (initially):**
- Large genome mapping optimizations
- RABEMA benchmarking mode
- Full SeqAn library port

**Estimated Total Effort:** 12-16 weeks for production-ready implementation

---

## Project Structure

```
yara-rs/
├── Cargo.toml
├── src/
│   ├── lib.rs                 # Library root
│   ├── main.rs                # CLI entry point
│   │
│   ├── index/
│   │   ├── mod.rs
│   │   ├── fm_index.rs        # FM-Index implementation
│   │   ├── suffix_array.rs    # SA construction
│   │   ├── bwt.rs             # Burrows-Wheeler Transform
│   │   └── serialization.rs   # Index I/O
│   │
│   ├── alignment/
│   │   ├── mod.rs
│   │   ├── seed.rs            # Seed extraction
│   │   ├── extend.rs          # Seed extension (Myers)
│   │   ├── verify.rs          # Match verification
│   │   └── cigar.rs           # CIGAR string generation
│   │
│   ├── mapper/
│   │   ├── mod.rs
│   │   ├── config.rs          # Configuration structs
│   │   ├── pipeline.rs        # Main mapping pipeline
│   │   ├── collector.rs       # Seed collection
│   │   ├── classifier.rs      # Read classification
│   │   ├── ranker.rs          # Seed ranking
│   │   ├── filter.rs          # Seed filtering
│   │   ├── extender.rs        # Hit extension
│   │   ├── writer.rs          # Output generation
│   │   └── stats.rs           # Performance statistics
│   │
│   ├── io/
│   │   ├── mod.rs
│   │   ├── fasta.rs           # FASTA/FASTQ parsing
│   │   ├── sam.rs             # SAM output
│   │   ├── bam.rs             # BAM output (via noodles)
│   │   └── prefetch.rs        # Async I/O
│   │
│   ├── types/
│   │   ├── mod.rs
│   │   ├── dna.rs             # DNA alphabet
│   │   ├── quality.rs         # Quality scores
│   │   ├── match_record.rs    # Match representation
│   │   └── hit.rs             # Hit representation
│   │
│   └── utils/
│       ├── mod.rs
│       ├── bits.rs            # Bit manipulation
│       └── parallel.rs        # Threading utilities
│
├── benches/
│   └── mapping_benchmark.rs
│
└── tests/
    ├── integration/
    │   └── hla_mapping.rs
    └── golden/
        └── (test data)
```

---

## Phase 1: Foundation (Weeks 1-3)

### 1.1 Project Setup and Dependencies

**Cargo.toml:**
```toml
[package]
name = "yara-rs"
version = "0.1.0"
edition = "2021"
license = "BSD-3-Clause"
description = "Fast DNA read aligner for HLA typing"

[dependencies]
# I/O
noodles = { version = "0.68", features = ["sam", "bam", "fasta", "fastq", "bgzf"] }
flate2 = "1.0"
bzip2 = "0.4"

# Parallelism
rayon = "1.8"
crossbeam = "0.8"

# Data structures
bio = "1.5"              # Bioinformatics primitives
packed_simd_2 = "0.3"    # SIMD operations (optional)

# CLI
clap = { version = "4.4", features = ["derive"] }

# Utilities
thiserror = "1.0"
anyhow = "1.0"
log = "0.4"
env_logger = "0.10"
memmap2 = "0.9"          # Memory-mapped I/O

[dev-dependencies]
criterion = "0.5"
proptest = "1.4"

[profile.release]
lto = true
codegen-units = 1
panic = "abort"

[[bin]]
name = "yara_indexer"
path = "src/bin/indexer.rs"

[[bin]]
name = "yara_mapper"
path = "src/bin/mapper.rs"
```

### 1.2 DNA Alphabet and Quality Types

**src/types/dna.rs:**
```rust
/// 2-bit encoded DNA base (A=0, C=1, G=2, T=3)
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Dna4 {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

/// 3-bit encoded DNA with N (A=0, C=1, G=2, T=3, N=4)
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Dna5 {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4,
}

impl Dna5 {
    #[inline]
    pub fn from_ascii(c: u8) -> Self {
        match c {
            b'A' | b'a' => Dna5::A,
            b'C' | b'c' => Dna5::C,
            b'G' | b'g' => Dna5::G,
            b'T' | b't' => Dna5::T,
            _ => Dna5::N,
        }
    }

    #[inline]
    pub fn complement(self) -> Self {
        match self {
            Dna5::A => Dna5::T,
            Dna5::T => Dna5::A,
            Dna5::C => Dna5::G,
            Dna5::G => Dna5::C,
            Dna5::N => Dna5::N,
        }
    }
}

/// Packed DNA sequence (4 bases per byte)
pub struct PackedDnaSeq {
    data: Vec<u8>,
    len: usize,
}
```

### 1.3 Match and Hit Types

**src/types/match_record.rs:**
```rust
/// Compact match representation (mirrors Yara's bit-packed Match)
#[derive(Clone, Copy, Debug)]
pub struct Match {
    /// Packed fields for memory efficiency
    /// Layout: [read_id: 24][contig_id: 16][contig_begin: 32][contig_end: 32][errors: 8][flags: 8]
    data: u128,
}

impl Match {
    pub fn new(
        read_id: u32,
        contig_id: u16,
        contig_begin: u32,
        contig_end: u32,
        errors: u8,
        reverse_strand: bool,
    ) -> Self {
        let flags = if reverse_strand { 1u8 } else { 0u8 };
        let data = (read_id as u128)
            | ((contig_id as u128) << 24)
            | ((contig_begin as u128) << 40)
            | ((contig_end as u128) << 72)
            | ((errors as u128) << 104)
            | ((flags as u128) << 112);
        Self { data }
    }

    #[inline]
    pub fn read_id(&self) -> u32 {
        (self.data & 0xFFFFFF) as u32
    }

    #[inline]
    pub fn contig_id(&self) -> u16 {
        ((self.data >> 24) & 0xFFFF) as u16
    }

    #[inline]
    pub fn contig_begin(&self) -> u32 {
        ((self.data >> 40) & 0xFFFFFFFF) as u32
    }

    #[inline]
    pub fn contig_end(&self) -> u32 {
        ((self.data >> 72) & 0xFFFFFFFF) as u32
    }

    #[inline]
    pub fn errors(&self) -> u8 {
        ((self.data >> 104) & 0xFF) as u8
    }

    #[inline]
    pub fn is_reverse_strand(&self) -> bool {
        ((self.data >> 112) & 1) != 0
    }

    /// Sort key for coordinate-based sorting
    #[inline]
    pub fn sort_key_coord(&self) -> u64 {
        ((self.contig_id() as u64) << 48)
            | ((self.is_reverse_strand() as u64) << 47)
            | ((self.contig_begin() as u64) << 15)
            | (self.errors() as u64)
    }
}
```

---

## Phase 2: FM-Index Implementation (Weeks 3-5)

### 2.1 FM-Index Core

**src/index/fm_index.rs:**
```rust
use memmap2::Mmap;
use std::fs::File;
use std::path::Path;

/// FM-Index for DNA sequences
pub struct FmIndex {
    /// BWT of the concatenated reference
    bwt: Vec<Dna5>,
    /// Occurrence table (sampled)
    occ: OccurrenceTable,
    /// Suffix array (sampled)
    sa: SampledSuffixArray,
    /// C table (cumulative character counts)
    c_table: [u64; 5],
    /// Contig boundaries in concatenated sequence
    contig_bounds: Vec<u64>,
    /// Contig names
    contig_names: Vec<String>,
}

impl FmIndex {
    /// Build index from FASTA file
    pub fn build<P: AsRef<Path>>(fasta_path: P) -> Result<Self, IndexError> {
        // 1. Read and concatenate sequences
        let (concat_seq, contig_bounds, contig_names) = Self::read_fasta(fasta_path)?;

        // 2. Build suffix array
        let sa = SuffixArray::new(&concat_seq);

        // 3. Build BWT from SA
        let bwt = Self::build_bwt(&concat_seq, &sa);

        // 4. Build occurrence table
        let occ = OccurrenceTable::new(&bwt);

        // 5. Build C table
        let c_table = Self::build_c_table(&bwt);

        // 6. Sample suffix array
        let sampled_sa = SampledSuffixArray::new(sa, 16); // Sample every 16

        Ok(Self {
            bwt,
            occ,
            sa: sampled_sa,
            c_table,
            contig_bounds,
            contig_names,
        })
    }

    /// Load index from disk (memory-mapped)
    pub fn load<P: AsRef<Path>>(index_prefix: P) -> Result<Self, IndexError> {
        // Memory-map index files for fast loading
        todo!()
    }

    /// Save index to disk
    pub fn save<P: AsRef<Path>>(&self, index_prefix: P) -> Result<(), IndexError> {
        todo!()
    }

    /// Backward search for exact pattern match
    /// Returns (start, end) range in suffix array
    pub fn backward_search(&self, pattern: &[Dna5]) -> Option<(u64, u64)> {
        if pattern.is_empty() {
            return None;
        }

        let mut start = 0u64;
        let mut end = self.bwt.len() as u64;

        for &c in pattern.iter().rev() {
            let c_idx = c as usize;
            start = self.c_table[c_idx] + self.occ.rank(c, start);
            end = self.c_table[c_idx] + self.occ.rank(c, end);

            if start >= end {
                return None;
            }
        }

        Some((start, end))
    }

    /// Get suffix array value at position (handles sampling)
    pub fn sa_value(&self, pos: u64) -> u64 {
        self.sa.get(pos, &self.bwt, &self.occ, &self.c_table)
    }

    /// Convert linear position to (contig_id, offset)
    pub fn to_contig_pos(&self, linear_pos: u64) -> (u16, u32) {
        let contig_id = self.contig_bounds
            .binary_search(&linear_pos)
            .unwrap_or_else(|i| i - 1);
        let offset = linear_pos - self.contig_bounds[contig_id];
        (contig_id as u16, offset as u32)
    }
}

/// Sampled occurrence table for space efficiency
pub struct OccurrenceTable {
    /// Sampled counts (every SAMPLE_RATE positions)
    samples: Vec<[u64; 5]>,
    /// Raw BWT for computing between samples
    bwt_blocks: Vec<u64>, // Packed 32 bases per u64
    sample_rate: usize,
}

impl OccurrenceTable {
    pub fn new(bwt: &[Dna5]) -> Self {
        const SAMPLE_RATE: usize = 64;
        // Implementation...
        todo!()
    }

    /// Count occurrences of character c in BWT[0..pos)
    #[inline]
    pub fn rank(&self, c: Dna5, pos: u64) -> u64 {
        // Get sampled count
        let sample_idx = (pos as usize) / self.sample_rate;
        let mut count = self.samples[sample_idx][c as usize];

        // Add counts from blocks between sample and pos
        // Use popcount for packed representation
        todo!()
    }
}
```

### 2.2 Approximate Search with Backtracking

**src/index/fm_index.rs (continued):**
```rust
impl FmIndex {
    /// Search with up to k errors using backtracking
    pub fn search_approximate<F>(
        &self,
        pattern: &[Dna5],
        max_errors: u8,
        mut callback: F,
    ) where
        F: FnMut(u64, u64, u8), // (sa_start, sa_end, errors)
    {
        let mut stack = Vec::with_capacity(pattern.len() * 4);

        // Initial state: full SA range, start from end of pattern
        stack.push(SearchState {
            sa_range: (0, self.bwt.len() as u64),
            pattern_pos: pattern.len(),
            errors: 0,
        });

        while let Some(state) = stack.pop() {
            if state.pattern_pos == 0 {
                // Found a match
                callback(state.sa_range.0, state.sa_range.1, state.errors);
                continue;
            }

            let pattern_char = pattern[state.pattern_pos - 1];

            // Try all characters
            for c in [Dna5::A, Dna5::C, Dna5::G, Dna5::T] {
                let new_start = self.c_table[c as usize]
                    + self.occ.rank(c, state.sa_range.0);
                let new_end = self.c_table[c as usize]
                    + self.occ.rank(c, state.sa_range.1);

                if new_start >= new_end {
                    continue;
                }

                let is_match = c == pattern_char || pattern_char == Dna5::N;
                let new_errors = state.errors + if is_match { 0 } else { 1 };

                if new_errors <= max_errors {
                    stack.push(SearchState {
                        sa_range: (new_start, new_end),
                        pattern_pos: state.pattern_pos - 1,
                        errors: new_errors,
                    });
                }
            }
        }
    }
}

struct SearchState {
    sa_range: (u64, u64),
    pattern_pos: usize,
    errors: u8,
}
```

---

## Phase 3: Alignment Algorithms (Weeks 5-7)

### 3.1 Myers Bit-Vector Algorithm

**src/alignment/extend.rs:**
```rust
/// Myers bit-vector algorithm for semi-global alignment
pub struct MyersAligner {
    /// Pattern length
    m: usize,
    /// Preprocessed pattern masks (one per character)
    peq: [u64; 5],
}

impl MyersAligner {
    pub fn new(pattern: &[Dna5]) -> Self {
        assert!(pattern.len() <= 64, "Pattern too long for single-word Myers");

        let m = pattern.len();
        let mut peq = [0u64; 5];

        for (i, &c) in pattern.iter().enumerate() {
            peq[c as usize] |= 1u64 << i;
        }

        Self { m, peq }
    }

    /// Find all alignments of pattern in text with <= max_errors
    /// Returns Vec<(end_position, errors)>
    pub fn find_all(&self, text: &[Dna5], max_errors: u8) -> Vec<(usize, u8)> {
        let mut results = Vec::new();

        // Initialize bit vectors
        let mut pv = !0u64; // All 1s
        let mut mv = 0u64;  // All 0s
        let mut score = self.m as i32;

        let high_bit = 1u64 << (self.m - 1);

        for (j, &c) in text.iter().enumerate() {
            let eq = self.peq[c as usize];

            let xv = eq | mv;
            let xh = ((eq & pv).wrapping_add(pv)) ^ pv | eq;

            let ph = mv | !(xh | pv);
            let mh = pv & xh;

            // Update score
            if (ph & high_bit) != 0 {
                score += 1;
            }
            if (mh & high_bit) != 0 {
                score -= 1;
            }

            // Shift for next column
            let ph_shift = (ph << 1) | 1;
            let mh_shift = mh << 1;

            pv = mh_shift | !(xv | ph_shift);
            mv = ph_shift & xv;

            if score <= max_errors as i32 {
                results.push((j + 1, score as u8));
            }
        }

        results
    }
}

/// Banded Myers for longer patterns
pub struct BandedMyersAligner {
    // For patterns > 64 bases, use word-parallel version
    words: usize,
    peq: Vec<[u64; 5]>,
}
```

### 3.2 Seed Extension

**src/alignment/extend.rs (continued):**
```rust
/// Seed extender using LCP optimization + Myers
pub struct SeedExtender<'a> {
    reference: &'a [Vec<Dna5>],  // Contigs
}

impl<'a> SeedExtender<'a> {
    pub fn new(reference: &'a [Vec<Dna5>]) -> Self {
        Self { reference }
    }

    /// Extend seed hit to full alignment
    pub fn extend(
        &self,
        read: &[Dna5],
        contig_id: u16,
        seed_contig_begin: u32,
        seed_contig_end: u32,
        seed_read_begin: usize,
        seed_read_end: usize,
        seed_errors: u8,
        max_errors: u8,
    ) -> Option<Match> {
        let contig = &self.reference[contig_id as usize];
        let mut total_errors = seed_errors;

        // LCP extension left
        let (left_contig_pos, left_errors) = self.extend_left(
            read,
            contig,
            seed_read_begin,
            seed_contig_begin as usize,
            max_errors - total_errors,
        )?;
        total_errors += left_errors;

        // LCP extension right
        let (right_contig_pos, right_errors) = self.extend_right(
            read,
            contig,
            seed_read_end,
            seed_contig_end as usize,
            max_errors - total_errors,
        )?;
        total_errors += right_errors;

        if total_errors <= max_errors {
            Some(Match::new(
                0, // read_id filled by caller
                contig_id,
                left_contig_pos as u32,
                right_contig_pos as u32,
                total_errors,
                false, // strand filled by caller
            ))
        } else {
            None
        }
    }

    /// Extend left using LCP then Myers
    fn extend_left(
        &self,
        read: &[Dna5],
        contig: &[Dna5],
        read_pos: usize,
        contig_pos: usize,
        max_errors: u8,
    ) -> Option<(usize, u8)> {
        if read_pos == 0 {
            return Some((contig_pos, 0));
        }

        // LCP trick: find exact match prefix
        let read_left = &read[..read_pos];
        let contig_start = contig_pos.saturating_sub(read_pos + max_errors as usize);
        let contig_left = &contig[contig_start..contig_pos];

        // Compute LCP from right
        let lcp = Self::lcp_reverse(read_left, contig_left);

        if lcp == read_pos {
            // Perfect match
            return Some((contig_pos - lcp, 0));
        }

        if max_errors == 0 {
            return None;
        }

        // Use Myers for remaining portion
        let read_remaining = &read[..read_pos - lcp];
        let contig_remaining = &contig[contig_start..contig_pos - lcp];

        // Reverse for left extension
        let read_rev: Vec<_> = read_remaining.iter().rev().copied().collect();
        let contig_rev: Vec<_> = contig_remaining.iter().rev().copied().collect();

        let aligner = MyersAligner::new(&read_rev);
        let hits = aligner.find_all(&contig_rev, max_errors);

        // Find best hit
        hits.into_iter()
            .filter(|(_, e)| *e <= max_errors)
            .min_by_key(|(_, e)| *e)
            .map(|(pos, errors)| (contig_pos - lcp - pos, errors))
    }

    fn extend_right(
        &self,
        read: &[Dna5],
        contig: &[Dna5],
        read_pos: usize,
        contig_pos: usize,
        max_errors: u8,
    ) -> Option<(usize, u8)> {
        // Similar to extend_left but forward direction
        todo!()
    }

    #[inline]
    fn lcp_reverse(a: &[Dna5], b: &[Dna5]) -> usize {
        a.iter()
            .rev()
            .zip(b.iter().rev())
            .take_while(|(x, y)| x == y)
            .count()
    }
}
```

---

## Phase 4: Mapping Pipeline (Weeks 7-10)

### 4.1 Pipeline Configuration

**src/mapper/config.rs:**
```rust
use clap::Parser;

#[derive(Parser, Debug, Clone)]
#[command(name = "yara_mapper", about = "Fast DNA read aligner")]
pub struct MapperConfig {
    /// Reference index prefix
    #[arg(required = true)]
    pub index: String,

    /// Input reads (FASTQ/FASTA, gzipped OK)
    #[arg(required = true)]
    pub reads: Vec<String>,

    /// Output file (SAM/BAM)
    #[arg(short, long, default_value = "-")]
    pub output: String,

    /// Error rate (0.0-1.0)
    #[arg(short, long, default_value = "0.05")]
    pub error_rate: f32,

    /// Number of threads
    #[arg(short, long, default_value = "0")]
    pub threads: usize,

    /// Secondary alignments mode: tag, record, omit
    #[arg(long, default_value = "tag")]
    pub secondary: SecondaryMode,

    /// Sensitivity: low, high, full
    #[arg(long, default_value = "high")]
    pub sensitivity: Sensitivity,

    /// Read group ID
    #[arg(long, default_value = "none")]
    pub read_group: String,

    /// Verbosity level
    #[arg(short, long, action = clap::ArgAction::Count)]
    pub verbose: u8,
}

#[derive(Clone, Copy, Debug, Default, clap::ValueEnum)]
pub enum SecondaryMode {
    #[default]
    Tag,
    Record,
    Omit,
}

#[derive(Clone, Copy, Debug, Default, clap::ValueEnum)]
pub enum Sensitivity {
    Low,
    #[default]
    High,
    Full,
}

impl MapperConfig {
    pub fn max_errors(&self, read_len: usize) -> u8 {
        ((read_len as f32) * self.error_rate).floor() as u8
    }

    pub fn thread_count(&self) -> usize {
        if self.threads == 0 {
            rayon::current_num_threads()
        } else {
            self.threads
        }
    }
}
```

### 4.2 Main Pipeline

**src/mapper/pipeline.rs:**
```rust
use rayon::prelude::*;
use crossbeam::channel;
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct Mapper {
    config: MapperConfig,
    index: FmIndex,
    stats: MapperStats,
}

impl Mapper {
    pub fn new(config: MapperConfig) -> Result<Self, MapperError> {
        let index = FmIndex::load(&config.index)?;
        Ok(Self {
            config,
            index,
            stats: MapperStats::default(),
        })
    }

    pub fn run(&mut self) -> Result<(), MapperError> {
        // Set up thread pool
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.config.thread_count())
            .build_global()?;

        // Open output
        let mut writer = SamWriter::new(&self.config.output, &self.index)?;
        writer.write_header(&self.config)?;

        // Process reads in batches
        let batch_size = 100_000;
        let mut reader = ReadBatchReader::new(&self.config.reads)?;

        while let Some(batch) = reader.next_batch(batch_size)? {
            let results = self.map_batch(&batch);
            writer.write_batch(&results)?;
        }

        writer.finish()?;

        if self.config.verbose > 0 {
            self.stats.print();
        }

        Ok(())
    }

    fn map_batch(&self, reads: &ReadBatch) -> Vec<MappingResult> {
        reads
            .par_iter()
            .map(|read| self.map_read(read))
            .collect()
    }

    fn map_read(&self, read: &Read) -> MappingResult {
        let max_errors = self.config.max_errors(read.seq.len());

        // 1. Extract seeds
        let seeds = self.extract_seeds(read, max_errors);

        // 2. Collect hits for each seed
        let hits: Vec<Hit> = seeds
            .iter()
            .flat_map(|seed| self.collect_hits(seed, max_errors))
            .collect();

        // 3. Rank and filter seeds/hits
        let ranked_hits = self.rank_hits(&hits, &seeds);

        // 4. Extend hits to matches
        let mut matches: Vec<Match> = ranked_hits
            .iter()
            .filter_map(|hit| self.extend_hit(read, hit, max_errors))
            .collect();

        // 5. Remove duplicates and sort
        self.deduplicate_matches(&mut matches);

        // 6. Compute mapping quality
        let mapq = self.compute_mapq(&matches, read.seq.len());

        MappingResult {
            read_id: read.id,
            matches,
            mapq,
        }
    }

    fn extract_seeds(&self, read: &Read, max_errors: u8) -> Vec<Seed> {
        let seed_len = Self::compute_seed_length(read.seq.len(), max_errors);
        let num_seeds = (read.seq.len() / seed_len) + 1;

        (0..num_seeds)
            .map(|i| {
                let start = i * seed_len;
                let end = (start + seed_len).min(read.seq.len());
                Seed {
                    read_begin: start,
                    read_end: end,
                    seq: read.seq[start..end].to_vec(),
                }
            })
            .collect()
    }

    fn compute_seed_length(read_len: usize, max_errors: u8) -> usize {
        // q-gram lemma: seed_len = read_len / (errors + 1)
        read_len / (max_errors as usize + 1)
    }

    fn collect_hits(&self, seed: &Seed, seed_errors: u8) -> Vec<Hit> {
        let mut hits = Vec::new();

        self.index.search_approximate(&seed.seq, seed_errors, |sa_start, sa_end, errors| {
            hits.push(Hit {
                sa_range: (sa_start, sa_end),
                seed_read_begin: seed.read_begin,
                seed_read_end: seed.read_end,
                errors,
            });
        });

        hits
    }

    fn extend_hit(&self, read: &Read, hit: &Hit, max_errors: u8) -> Option<Match> {
        // For each SA position in the hit range
        for sa_pos in hit.sa_range.0..hit.sa_range.1 {
            let linear_pos = self.index.sa_value(sa_pos);
            let (contig_id, contig_offset) = self.index.to_contig_pos(linear_pos);

            // Try extension
            let extender = SeedExtender::new(&self.index.contigs());
            if let Some(mut m) = extender.extend(
                &read.seq,
                contig_id,
                contig_offset,
                contig_offset + (hit.seed_read_end - hit.seed_read_begin) as u32,
                hit.seed_read_begin,
                hit.seed_read_end,
                hit.errors,
                max_errors,
            ) {
                // Fill in read ID
                m.set_read_id(read.id);
                return Some(m);
            }
        }
        None
    }

    fn deduplicate_matches(&self, matches: &mut Vec<Match>) {
        matches.sort_by_key(|m| m.sort_key_coord());
        matches.dedup_by(|a, b| {
            a.contig_id() == b.contig_id()
                && a.contig_begin() == b.contig_begin()
                && a.is_reverse_strand() == b.is_reverse_strand()
        });
    }

    fn compute_mapq(&self, matches: &[Match], read_len: usize) -> u8 {
        if matches.is_empty() {
            return 0;
        }

        let best_errors = matches.iter().map(|m| m.errors()).min().unwrap();
        let best_count = matches.iter().filter(|m| m.errors() == best_errors).count();
        let sub_count = matches.len() - best_count;

        let error_rate = best_errors as f64 / read_len as f64;
        let prob = Self::match_probability(error_rate, best_count, sub_count);

        (-10.0 * (1.0 - prob.min(0.9999999)).log10()).round() as u8
    }

    fn match_probability(error_rate: f64, best_count: usize, sub_count: usize) -> f64 {
        // Simplified MAPQ calculation
        // Full implementation would use proper probabilistic model
        let p_correct = 1.0 / (best_count as f64 + 0.1 * sub_count as f64);
        p_correct * (1.0 - error_rate).powi(100)
    }
}
```

### 4.3 SAM/BAM Output

**src/io/sam.rs:**
```rust
use noodles::sam::{self, Header, Record};
use std::io::Write;

pub struct SamWriter {
    inner: Box<dyn Write>,
    header: Header,
    format: OutputFormat,
}

impl SamWriter {
    pub fn new(path: &str, index: &FmIndex) -> Result<Self, IoError> {
        let (inner, format): (Box<dyn Write>, OutputFormat) = if path == "-" {
            (Box::new(std::io::stdout()), OutputFormat::Sam)
        } else if path.ends_with(".bam") {
            let file = std::fs::File::create(path)?;
            // Use noodles BAM writer
            todo!()
        } else {
            let file = std::fs::File::create(path)?;
            (Box::new(std::io::BufWriter::new(file)), OutputFormat::Sam)
        };

        let header = Self::build_header(index);

        Ok(Self { inner, header, format })
    }

    pub fn write_header(&mut self, config: &MapperConfig) -> Result<(), IoError> {
        // Write SAM header
        writeln!(self.inner, "@HD\tVN:1.4\tSO:unsorted")?;

        for (id, name) in self.header.reference_sequences().iter().enumerate() {
            writeln!(self.inner, "@SQ\tSN:{}\tLN:{}", name.0, name.1)?;
        }

        writeln!(self.inner, "@PG\tID:yara-rs\tPN:yara-rs\tVN:{}", env!("CARGO_PKG_VERSION"))?;
        writeln!(self.inner, "@RG\tID:{}\tSM:{}", config.read_group, config.read_group)?;

        Ok(())
    }

    pub fn write_result(&mut self, result: &MappingResult, read: &Read) -> Result<(), IoError> {
        if result.matches.is_empty() {
            self.write_unmapped(read)?;
        } else {
            self.write_mapped(result, read)?;
        }
        Ok(())
    }

    fn write_mapped(&mut self, result: &MappingResult, read: &Read) -> Result<(), IoError> {
        let primary = &result.matches[0];

        // Primary alignment
        let flag = if primary.is_reverse_strand() { 0x10 } else { 0 };

        write!(
            self.inner,
            "{}\t{}\t{}\t{}\t{}\t{}M\t*\t0\t0\t{}\t{}\tNM:i:{}\tX0:i:{}\tX1:i:{}",
            read.name,
            flag,
            self.header.reference_sequences().get_index(primary.contig_id() as usize).unwrap().0,
            primary.contig_begin() + 1,
            result.mapq,
            read.seq.len(),
            Self::seq_to_string(&read.seq, primary.is_reverse_strand()),
            Self::qual_to_string(&read.qual),
            primary.errors(),
            result.matches.iter().filter(|m| m.errors() == primary.errors()).count(),
            result.matches.len() - result.matches.iter().filter(|m| m.errors() == primary.errors()).count(),
        )?;

        // XA tag for secondary alignments
        if result.matches.len() > 1 {
            write!(self.inner, "\tXA:Z:")?;
            for m in result.matches.iter().skip(1) {
                write!(
                    self.inner,
                    "{},{},{},{},{};",
                    self.header.reference_sequences().get_index(m.contig_id() as usize).unwrap().0,
                    m.contig_begin() + 1,
                    m.contig_end() + 1,
                    if m.is_reverse_strand() { '-' } else { '+' },
                    m.errors(),
                )?;
            }
        }

        writeln!(self.inner)?;
        Ok(())
    }

    fn write_unmapped(&mut self, read: &Read) -> Result<(), IoError> {
        writeln!(
            self.inner,
            "{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}",
            read.name,
            Self::seq_to_string(&read.seq, false),
            Self::qual_to_string(&read.qual),
        )
    }

    fn seq_to_string(seq: &[Dna5], reverse: bool) -> String {
        if reverse {
            seq.iter().rev().map(|b| b.complement().to_ascii()).collect()
        } else {
            seq.iter().map(|b| b.to_ascii()).collect()
        }
    }

    fn qual_to_string(qual: &[u8]) -> String {
        qual.iter().map(|&q| (q + 33) as char).collect()
    }
}
```

---

## Phase 5: Testing and Validation (Weeks 10-12)

### 5.1 Unit Tests

**src/index/fm_index.rs:**
```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_backward_search_exact() {
        let seq = "ACGTACGT".chars().map(|c| Dna5::from_ascii(c as u8)).collect::<Vec<_>>();
        let index = FmIndex::build_from_seq(&seq);

        // Search for "ACG"
        let pattern: Vec<Dna5> = "ACG".chars().map(|c| Dna5::from_ascii(c as u8)).collect();
        let result = index.backward_search(&pattern);

        assert!(result.is_some());
        let (start, end) = result.unwrap();
        assert_eq!(end - start, 2); // Two occurrences
    }

    #[test]
    fn test_approximate_search() {
        let seq = "ACGTACGT".chars().map(|c| Dna5::from_ascii(c as u8)).collect::<Vec<_>>();
        let index = FmIndex::build_from_seq(&seq);

        let pattern: Vec<Dna5> = "ACT".chars().map(|c| Dna5::from_ascii(c as u8)).collect();
        let mut hits = Vec::new();
        index.search_approximate(&pattern, 1, |s, e, err| hits.push((s, e, err)));

        // Should find "ACG" with 1 error
        assert!(!hits.is_empty());
    }
}
```

### 5.2 Integration Tests

**tests/integration/hla_mapping.rs:**
```rust
use yara_rs::*;
use std::path::PathBuf;

#[test]
fn test_hla_reference_indexing() {
    let test_data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let hla_fasta = test_data.join("hla_ref.fa");

    let index = FmIndex::build(&hla_fasta).expect("Failed to build index");

    // Save and reload
    let index_path = test_data.join("hla_ref.idx");
    index.save(&index_path).expect("Failed to save index");

    let loaded = FmIndex::load(&index_path).expect("Failed to load index");
    assert_eq!(index.contig_count(), loaded.contig_count());
}

#[test]
fn test_mapping_output_matches_original() {
    // Compare output with original Yara C++ implementation
    let test_data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");

    let config = MapperConfig {
        index: test_data.join("hla_ref").to_string_lossy().to_string(),
        reads: vec![test_data.join("reads.fq").to_string_lossy().to_string()],
        output: test_data.join("output.sam").to_string_lossy().to_string(),
        ..Default::default()
    };

    let mut mapper = Mapper::new(config).expect("Failed to create mapper");
    mapper.run().expect("Mapping failed");

    // Compare SAM output (excluding @PG line)
    let expected = std::fs::read_to_string(test_data.join("expected.sam")).unwrap();
    let actual = std::fs::read_to_string(test_data.join("output.sam")).unwrap();

    let expected_records: Vec<_> = expected.lines().filter(|l| !l.starts_with("@PG")).collect();
    let actual_records: Vec<_> = actual.lines().filter(|l| !l.starts_with("@PG")).collect();

    assert_eq!(expected_records, actual_records);
}
```

### 5.3 Benchmarks

**benches/mapping_benchmark.rs:**
```rust
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use yara_rs::*;

fn bench_indexing(c: &mut Criterion) {
    let test_data = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");

    c.bench_function("index_hla_reference", |b| {
        b.iter(|| {
            FmIndex::build(test_data.join("hla_ref.fa")).unwrap()
        })
    });
}

fn bench_mapping(c: &mut Criterion) {
    let test_data = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let index = FmIndex::load(test_data.join("hla_ref")).unwrap();

    let mut group = c.benchmark_group("mapping_throughput");

    for thread_count in [1, 2, 4, 8] {
        group.bench_with_input(
            BenchmarkId::new("threads", thread_count),
            &thread_count,
            |b, &threads| {
                let config = MapperConfig {
                    threads,
                    ..Default::default()
                };
                b.iter(|| {
                    // Map reads
                })
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_indexing, bench_mapping);
criterion_main!(benches);
```

---

## Phase 6: OptiType Integration (Weeks 12-14)

### 6.1 OptiType-Specific Optimizations

**src/mapper/optitype.rs:**
```rust
/// OptiType-optimized mapping configuration
pub fn optitype_config() -> MapperConfig {
    MapperConfig {
        error_rate: 0.05,
        secondary: SecondaryMode::Tag,  // OptiType uses XA tags
        sensitivity: Sensitivity::High,
        ..Default::default()
    }
}

/// Pre-built HLA reference index
pub struct HlaIndex {
    index: FmIndex,
    // HLA-specific metadata
    allele_groups: HashMap<String, Vec<u16>>,  // HLA-A, HLA-B, etc.
}

impl HlaIndex {
    /// Load IMGT/HLA reference
    pub fn load_imgt(version: &str) -> Result<Self, IndexError> {
        todo!()
    }

    /// Filter matches to specific HLA locus
    pub fn filter_locus(&self, matches: &[Match], locus: &str) -> Vec<Match> {
        let contig_ids = self.allele_groups.get(locus).unwrap_or(&vec![]);
        matches
            .iter()
            .filter(|m| contig_ids.contains(&m.contig_id()))
            .cloned()
            .collect()
    }
}
```

### 6.2 Output Compatibility

Ensure SAM output format matches what OptiType expects:
- XA tag format for alternative alignments
- Proper MAPQ calculation
- NM, X0, X1 tags

---

## Phase 7: Performance Optimization (Weeks 14-16)

### 7.1 SIMD Optimizations

**src/alignment/simd.rs:**
```rust
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// SIMD-accelerated LCP computation
#[cfg(target_arch = "x86_64")]
pub fn lcp_length_simd(a: &[u8], b: &[u8]) -> usize {
    let len = a.len().min(b.len());
    let mut i = 0;

    // Process 32 bytes at a time with AVX2
    #[cfg(target_feature = "avx2")]
    unsafe {
        while i + 32 <= len {
            let va = _mm256_loadu_si256(a[i..].as_ptr() as *const __m256i);
            let vb = _mm256_loadu_si256(b[i..].as_ptr() as *const __m256i);
            let cmp = _mm256_cmpeq_epi8(va, vb);
            let mask = _mm256_movemask_epi8(cmp) as u32;

            if mask != 0xFFFFFFFF {
                return i + mask.trailing_ones() as usize;
            }
            i += 32;
        }
    }

    // Scalar fallback
    while i < len && a[i] == b[i] {
        i += 1;
    }

    i
}
```

### 7.2 Memory Pool Allocator

**src/utils/pool.rs:**
```rust
use std::cell::RefCell;

thread_local! {
    static MATCH_POOL: RefCell<Vec<Match>> = RefCell::new(Vec::with_capacity(10000));
}

pub fn with_match_pool<F, R>(f: F) -> R
where
    F: FnOnce(&mut Vec<Match>) -> R,
{
    MATCH_POOL.with(|pool| {
        let mut pool = pool.borrow_mut();
        pool.clear();
        f(&mut pool)
    })
}
```

---

## Deliverables Checklist

### Week 3
- [ ] Project structure created
- [ ] DNA types implemented
- [ ] Match/Hit types implemented
- [ ] Basic FASTA reader

### Week 5
- [ ] FM-Index construction
- [ ] Index serialization
- [ ] Exact backward search

### Week 7
- [ ] Approximate search with backtracking
- [ ] Myers bit-vector alignment
- [ ] Seed extension

### Week 10
- [ ] Full mapping pipeline
- [ ] SAM output
- [ ] Multi-threading with rayon

### Week 12
- [ ] Unit test coverage > 80%
- [ ] Integration tests pass
- [ ] Output matches original Yara

### Week 14
- [ ] OptiType integration tested
- [ ] BAM output support

### Week 16
- [ ] SIMD optimizations
- [ ] Performance benchmarks
- [ ] Documentation complete
- [ ] Release build optimized

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| FM-Index correctness | Extensive test suite, compare with existing implementations |
| Performance regression | Benchmark against C++ version throughout development |
| Memory usage | Profile with large datasets, implement streaming where possible |
| Threading bugs | Use Rust's ownership system, minimize shared mutable state |
| OptiType compatibility | Test with real HLA typing workflows early |

---

## Dependencies on External Crates

| Crate | Purpose | Alternative |
|-------|---------|-------------|
| noodles | SAM/BAM I/O | rust-htslib (but larger dependency) |
| rayon | Parallelism | crossbeam + manual thread pool |
| bio | Bioinformatics utilities | Implement from scratch |
| memmap2 | Memory-mapped I/O | std::fs with manual buffering |
| clap | CLI parsing | structopt or manual |

---

## Success Criteria

1. **Functional parity:** Output matches original Yara for HLA calling
2. **Performance:** Within 20% of C++ version (ideally faster)
3. **Memory:** No more than 1.5x memory usage of C++ version
4. **Reliability:** No crashes on malformed input
5. **Usability:** Drop-in replacement for OptiType workflows
