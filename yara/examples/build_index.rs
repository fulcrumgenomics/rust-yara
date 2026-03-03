//! Example: build a YARA FM index from a FASTA reference file.
//!
//! Usage:
//!
//!     cargo run --example build_index -- <reference.fasta> [output-prefix]
//!
//! If the output prefix is omitted, the FASTA filename without extension is used.

use yara_mapper::{IndexerOptions, YaraIndexer};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <reference.fasta> [output_prefix]", args[0]);
        std::process::exit(1);
    }
    let fasta = &args[1];
    let prefix = if args.len() >= 3 {
        args[2].clone()
    } else {
        fasta.rsplit_once('.').map_or_else(|| fasta.clone(), |(base, _)| base.to_string())
    };

    let options = IndexerOptions { verbose: 1, ..Default::default() };

    eprintln!("Building index from: {fasta}");
    eprintln!("Output prefix: {prefix}");

    let indexer = YaraIndexer::build(fasta, &prefix, &options).unwrap_or_else(|e| {
        eprintln!("Failed to build index: {e}");
        std::process::exit(1);
    });

    eprintln!("Built index with {} contigs", indexer.contig_count());
    for (i, (name, len)) in indexer.contig_names().iter().zip(indexer.contig_lengths()).enumerate()
    {
        eprintln!("  contig {i}: {name} ({len} bp)");
        if i >= 5 {
            eprintln!("  ... ({} more)", indexer.contig_count() - i - 1);
            break;
        }
    }
}
