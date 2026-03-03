//! Safe Rust wrapper for the YARA read mapper.
//!
//! This crate provides a high-level API for building and loading YARA FM
//! indices, mapping paired-end reads in batches, and retrieving alignment
//! results as owned Rust types — without SAM serialization/deserialization
//! overhead.
//!
//! ## Example
//!
//! ```no_run
//! use yara_mapper::{IndexerOptions, MapperOptions, ReadBatch, ReadEnd, YaraIndexer, YaraMapper};
//!
//! // Build an index from a FASTA file.
//! let indexer = YaraIndexer::build("ref.fasta", "ref", &IndexerOptions::default()).unwrap();
//! println!("Indexed {} contigs", indexer.contig_count());
//!
//! // Open the index and map reads.
//! let mapper = YaraMapper::open("ref", &MapperOptions::default()).unwrap();
//! let mut batch = ReadBatch::new();
//! batch.push(
//!     "read1",
//!     ReadEnd { seq: b"ACGTACGT", qual: b"IIIIIIII" },
//!     ReadEnd { seq: b"TGCATGCA", qual: b"IIIIIIII" },
//! ).unwrap();
//! let records = mapper.map_paired(&batch).unwrap();
//! for rec in &records {
//!     println!("contig={} pos={} mapq={}", rec.contig_id, rec.pos, rec.mapq);
//! }
//! ```

mod error;
mod ffi_helpers;
mod indexer;
mod mapper;
mod options;
mod record;

pub use error::YaraError;
pub use indexer::{IndexerOptions, YaraIndexer};
pub use mapper::{ReadBatch, ReadEnd, YaraMapper};
pub use options::{MapperOptions, SecondaryMode, Sensitivity};
pub use record::{CigarOp, YaraRecord};
