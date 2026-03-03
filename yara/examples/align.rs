//! Example: align a small batch of paired-end reads using the YARA mapper.
//!
//! Usage:
//!
//!     cargo run --example align -- <index-prefix>
//!
//! The `index_prefix` should point to a pre-built YARA index (e.g.,
//! `ref/hla.fasta` if the index files are `ref/hla.fasta.yara.*`).

use yara_mapper::{MapperOptions, ReadBatch, ReadEnd, Sensitivity, YaraMapper};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <index_prefix>", args[0]);
        std::process::exit(1);
    }
    let index_prefix = &args[1];

    // Configure the mapper.
    let options = MapperOptions {
        sensitivity: Sensitivity::Full,
        threads: 2,
        verbose: 1,
        ..Default::default()
    };

    // Open the index.
    eprintln!("Opening index: {index_prefix}");
    let mapper = YaraMapper::open(index_prefix, &options).unwrap_or_else(|e| {
        eprintln!("Failed to open index: {e}");
        std::process::exit(1);
    });

    eprintln!("Loaded {} contigs", mapper.contig_count());
    for (i, (name, len)) in mapper.contig_names().iter().zip(mapper.contig_lengths()).enumerate() {
        eprintln!("  contig {i}: {name} ({len} bp)");
        if i >= 5 {
            eprintln!("  ... ({} more)", mapper.contig_count() - i - 1);
            break;
        }
    }

    // Create a synthetic read pair (these won't align to anything real).
    let mut batch = ReadBatch::new();
    batch
        .push(
            "test_read_1",
            ReadEnd {
                seq: b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                qual: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            },
            ReadEnd {
                seq: b"TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
                qual: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            },
        )
        .expect("failed to add read pair");

    // Map the reads.
    eprintln!("Mapping {} read pair(s)...", batch.len());
    match mapper.map_paired(&batch) {
        Ok(records) => {
            eprintln!("Got {} alignment record(s)", records.len());
            for rec in &records {
                println!(
                    "pair={} read1={} contig={} pos={} rev={} unmapped={} secondary={} mapq={} nm={} cigar={} flag={}",
                    rec.read_pair_index,
                    rec.is_read1,
                    rec.contig_id,
                    rec.pos,
                    rec.is_reverse,
                    rec.is_unmapped,
                    rec.is_secondary,
                    rec.mapq,
                    rec.nm,
                    rec.cigar_string(),
                    rec.flag,
                );
            }
        }
        Err(e) => {
            eprintln!("Mapping failed: {e}");
            std::process::exit(1);
        }
    }
}
