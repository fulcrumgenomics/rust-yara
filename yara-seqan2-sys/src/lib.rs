//! Vendored `SeqAn2` and YARA app headers.
//!
//! This crate contains no Rust code. It vendors the `SeqAn2` `include/` directory
//! and the YARA application headers (`apps/yara/`) and exposes their paths via
//! cargo `links` metadata (`DEP_SEQAN2_INCLUDE` and `DEP_SEQAN2_YARA_APP`).
