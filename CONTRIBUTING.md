# Contributing to rust-yara

## Development Setup

### Prerequisites

- Rust (stable toolchain — see `rust-toolchain.toml` for the pinned version)
- zlib development headers
- OpenMP:
  - macOS: `brew install libomp`
  - Linux: `libomp-dev` or `libgomp` (typically included with GCC)

### Running Checks Manually

```bash
# Format check (fails if formatting differs)
cargo ci-fmt

# Lint check (fails on any warnings)
cargo ci-lint

# Build all targets
cargo ci-build
```

> **Note:** rust-yara does not have tests yet.  Test infrastructure is planned
> for a future release.

## Code Style

- Run `cargo fmt` before committing
- Fix all clippy warnings
- Add backticks around identifiers in doc comments (e.g., `` `contig_count` ``)

## Pull Requests

1. Ensure all CI checks pass (`cargo ci-fmt`, `cargo ci-lint`, `cargo ci-build`)
2. Keep PRs focused and reasonably sized
3. Include tests for new functionality (when test infrastructure is available)
