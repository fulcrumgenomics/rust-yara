# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Yara (Yet Another Read Aligner) is an exact DNA sequence read aligner for mapping sequencing reads to reference genomes. It's part of the SeqAn library ecosystem.

**Two executables:**
- `yara_indexer` - Builds FM-index of reference genomes
- `yara_mapper` - Maps DNA reads to indexed reference genomes

## Build Commands

Build from the SeqAn root directory:
```bash
mkdir yara-build && cd yara-build
cmake ../seqan -DSEQAN_BUILD_SYSTEM=APP:yara
make all
```

**Dependencies:** OpenMP (threading), ZLIB (gzip/BAM), BZip2 (optional)

**Build option:** `-DYARA_LARGE_CONTIGS=ON` (default) enables support for >32k contigs or contigs >4Gbp

## Running Tests

```bash
python tests/run_tests.py SOURCE_ROOT_PATH BINARY_ROOT_PATH
```

The test framework uses SeqAn's `app_tests` Python module. Tests run the indexer on `tests/input/adeno-genome.fa` and mapper on `tests/input/adeno-reads_1.fa`, comparing outputs against golden files in `tests/gold/`.

## Architecture

### Mapper Pipeline (mapper_*.h files)

The mapper follows a staged pipeline pattern. Each stage is a separate header file:

```
SeedsCollector (mapper_collector.h)    - Extract candidate seeds from reads
       ↓
ReadsClassifier (mapper_classifier.h)  - Classify reads by seed count
       ↓
SeedsRanker (mapper_ranker.h)          - Rank seeds by quality
       ↓
SeedsFilter (mapper_filter.h)          - Filter low-quality seeds
       ↓
HitsExtender (mapper_extender.h)       - Extend seed hits to full matches
       ↓
MatchesVerifier (mapper_verifier.h)    - Verify match quality
       ↓
MatchesAligner (mapper_aligner.h)      - Generate alignments (CIGAR)
       ↓
MatchesWriter (mapper_writer.h)        - Output to SAM/BAM
```

### Key Files

- `mapper.cpp` / `mapper.h` - Main entry point and configuration
- `indexer.cpp` - Index builder entry point
- `bits_*.h` - Data structures (Match, Hit, ReadsContext) with bit-packed storage
- `find_extender.h` / `find_verifier.h` - Core alignment algorithms
- `index_fm.h` - FM-Index overloads

### Configuration System

Template-based configuration via `ReadMapperConfig`:
- `TThreading` - Parallel/Serial execution
- `TSequencing` - SingleEnd/PairedEnd mode
- `TSeedsDistance` - HammingDistance/EditDistance
- `TContigsSize/Len/Sum` - Type sizes for different reference genome scales

Configuration is resolved at compile-time through template chain in `mapper.cpp` (`configureMapper()` functions).

## Code Conventions

- Namespace: `seqan2`
- Header guards: `APP_YARA_[MODULE]_H_`
- Stage boundaries marked with `// ============` comment separators
- Heavy use of SeqAn library types (String, StringSet, Index)
