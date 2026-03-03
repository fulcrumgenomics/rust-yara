[![Build Status](https://github.com/fulcrumgenomics/rust-yara/actions/workflows/check.yml/badge.svg)](https://github.com/fulcrumgenomics/rust-yara/actions/workflows/check.yml)
[![crates.io](https://img.shields.io/crates/v/yara-mapper.svg)](https://crates.io/crates/yara-mapper)
[![docs.rs](https://docs.rs/yara-mapper/badge.svg)](https://docs.rs/yara-mapper)
[![License: MIT](https://img.shields.io/crates/l/yara-mapper.svg)](LICENSE)

# rust-yara

Rust FFI bindings for the YARA read mapper.

## Overview

rust-yara provides safe Rust bindings for the
[YARA read mapper](https://github.com/seqan/seqan/tree/main/apps/yara), a fast
and accurate read mapper based on the SeqAn2 C++ library.  The workspace
supports both building FM indices from FASTA reference files and mapping
paired-end reads against a pre-built index — all from Rust, without
SAM serialization/deserialization overhead.

## Crate Structure

| Crate | Description |
|-------|-------------|
| [`yara-seqan2-sys`](https://crates.io/crates/yara-seqan2-sys) | Vendored SeqAn2 headers for the YARA read mapper |
| [`yara-mapper-sys`](https://crates.io/crates/yara-mapper-sys) | Low-level FFI bindings to the YARA read mapper (SeqAn2) |
| [`yara-mapper`](https://crates.io/crates/yara-mapper) | Safe Rust wrapper for the YARA read mapper |

## Prerequisites

- **Rust** stable toolchain (see `rust-toolchain.toml` for the pinned version)
- **zlib** development headers
- **OpenMP**:
  - macOS: `brew install libomp`
  - Linux: `libgomp` (typically included with GCC; install `libomp-dev` if
    using Clang)

## Quick Start

Add `yara-mapper` to your `Cargo.toml`:

```sh
cargo add yara-mapper
```

Build an index and map reads:

```rust,no_run
use yara_mapper::{IndexerOptions, MapperOptions, ReadBatch, YaraIndexer, YaraMapper};

// Build an index from a FASTA file.
let indexer = YaraIndexer::build("ref.fasta", "ref", &IndexerOptions::default()).unwrap();
println!("Indexed {} contigs", indexer.contig_count());

// Open the index and map reads.
let mapper = YaraMapper::open("ref", &MapperOptions::default()).unwrap();
let mut batch = ReadBatch::new();
batch.push("read1", b"ACGTACGT", b"IIIIIIII", b"TGCATGCA", b"IIIIIIII").unwrap();
let records = mapper.map_paired(&batch).unwrap();
for rec in &records {
    println!("contig={} pos={} mapq={}", rec.contig_id, rec.pos, rec.mapq);
}
```

## Examples

Build an FM index from a FASTA reference:

```sh
cargo run --example build_index -- ref.fasta ref
```

Map reads against a pre-built index:

```sh
cargo run --example align -- ref
```

## Building from Source

```sh
cargo build                  # debug build
cargo clippy --all-targets   # lint
cargo fmt --check            # format check
```

To use a local SeqAn2 checkout instead of the vendored headers, set `SEQAN_DIR`:

```sh
export SEQAN_DIR=/path/to/seqan
```

## License

This project is licensed under the [MIT License](LICENSE).

The vendored [SeqAn2](https://github.com/seqan/seqan) headers and
[YARA](https://github.com/seqan/seqan/tree/main/apps/yara) application sources
are licensed under the BSD 3-Clause License — see
[`yara-seqan2-sys/vendor/include/seqan/LICENSE`](yara-seqan2-sys/vendor/include/seqan/LICENSE)
and
[`yara-seqan2-sys/vendor/apps/yara/LICENSE`](yara-seqan2-sys/vendor/apps/yara/LICENSE).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and guidelines.
