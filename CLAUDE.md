# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

SeqAn2 headers are vendored in the `yara-seqan2-sys` crate, so no external checkout is needed:

```bash
cargo build                         # debug build
cargo build --release               # release build
cargo clippy --all-targets          # lint (no clippy.toml; use default rules)
cargo fmt --check                   # format check
cargo run --example build_index -- <ref.fasta>  # build index (from yara-mapper crate)
cargo run --example align -- <index_prefix>   # map reads  (from yara-mapper crate)
```

To use a local SeqAn2 checkout instead of the vendored headers, set `SEQAN_DIR`:

```bash
export SEQAN_DIR=/path/to/seqan    # optional override (contains include/ and apps/yara/)
```

On macOS, OpenMP is detected via Homebrew libomp (`/opt/homebrew/opt/libomp/`). Override with `HOMEBREW_PREFIX` if needed. On Linux, libgomp is used.

## Architecture

Three-crate workspace wrapping the YARA read mapper and indexer (SeqAn2 C++ template library) via FFI, allowing Rust code to build an FM index from a FASTA file and map paired-end DNA reads without SAM serialization.

```
Rust caller
  → yara-mapper crate (safe API: YaraIndexer, YaraMapper, ReadBatch, YaraRecord)
    → yara-mapper-sys crate (raw extern "C" FFI bindings + build.rs)
      → yara-seqan2-sys crate (vendored SeqAn2 headers, exposes paths via links metadata)
      → yara_indexer_shim.cpp (C++ shim wrapping SeqAn2 indexer templates)
      → yara_shim.cpp (C++ shim wrapping SeqAn2 mapper templates)
        → YARA indexer / mapper (SeqAn2 header-only library)
```

### yara-seqan2-sys (vendored headers)

- `vendor/include/` — full SeqAn2 header-only library
- `vendor/apps/yara/` — YARA application headers
- `build.rs` — exposes `DEP_SEQAN2_INCLUDE` and `DEP_SEQAN2_YARA_APP` via cargo `links` metadata
- No Rust code; header-only, no compilation

### yara-mapper-sys (low-level)

- `build.rs` — compiles `cpp/yara_shim.cpp` and `cpp/yara_indexer_shim.cpp` via the `cc` crate, linking SeqAn2 headers, zlib, and OpenMP
- `cpp/yara_shim.h` / `cpp/yara_shim.cpp` — mapper C API: 7 functions, 4 `repr(C)` structs; `RecordCollector` replaces YARA's SAM writer; factory chain dispatches by contig limits + sensitivity
- `cpp/yara_indexer_shim.h` / `cpp/yara_indexer_shim.cpp` — indexer C API: 5 functions, 2 structs; builds FM index from FASTA, saves 5 index files (`.txt`, `.rid`, `.txt.size`, `.sa`, `.lf`); factory chain dispatches by contig limits
- `src/lib.rs` — raw `unsafe extern "C"` declarations matching both C headers

### yara-mapper (safe wrapper)

- `indexer.rs` — `YaraIndexer` (owns a `NonNull<YaraIndexerHandle>`, `Send` but not `Sync`), `IndexerOptions`
- `mapper.rs` — `YaraMapper` (owns a `NonNull<YaraMapperHandle>`, `Send` but not `Sync` due to internal OpenMP), `ReadBatch` (collects `CString` pairs)
- `options.rs` — `MapperOptions` with `to_ffi()` conversion; enums `Sensitivity` and `SecondaryMode`
- `record.rs` — `YaraRecord` (fully owned, all heap data copied from C++ then freed), `CigarOp`
- `error.rs` — `YaraError` enum: `IndexOpen`, `IndexBuild`, `Mapping`, `InvalidInput`

### Memory ownership across the FFI boundary

C++ allocates heap memory for CIGAR arrays, seq/qual strings, and XA tags inside each `YaraAlignmentRecord`. The Rust side copies these into owned types (`Vec<CigarOp>`, `Vec<u8>`, `String`) in `convert_record()`, then calls `yara_mapper_free_records()` to deallocate the C++ memory. On error (`map_paired` returns -1), C++ frees partial records itself before returning — the Rust side does not call `free_records` in the error path.

The output buffer is zero-initialized once by C++ at the start of `mapPaired` (single memset). Neither `_nextRecord` nor the Rust caller perform additional zeroing.

### Thread safety

`YaraMapper` is `Send` but not `Sync`. OpenMP parallelism is internal to each `map_paired` call; concurrent calls from multiple threads are not safe.
