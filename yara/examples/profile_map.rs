//! Profiling harness: load real FASTQ read pairs and map them in batches.
//!
//! Usage:
//!
//!     cargo run --release --example profile_map -- \
//!         <index-prefix> <r1.fastq.gz> <r2.fastq.gz> [batch_size] [max_reads]
//!
//! Designed to be run under `samply record` or `cargo flamegraph`.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

use flate2::read::MultiGzDecoder;
use yara_mapper::{MapperOptions, ReadBatch, ReadEnd, Sensitivity, YaraMapper};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!(
            "Usage: {} <index_prefix> <r1.fastq.gz> <r2.fastq.gz> [batch_size] [max_reads]",
            args[0]
        );
        std::process::exit(1);
    }
    let index_prefix = &args[1];
    let r1_path = &args[2];
    let r2_path = &args[3];
    let batch_size: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(1000);
    let max_reads: usize = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(usize::MAX);

    let options = MapperOptions {
        sensitivity: Sensitivity::Full,
        threads: 1,
        verbose: 0,
        ..Default::default()
    };

    eprintln!("Opening index: {index_prefix}");
    let mapper = YaraMapper::open(index_prefix, &options).unwrap_or_else(|e| {
        eprintln!("Failed to open index: {e}");
        std::process::exit(1);
    });
    eprintln!("Loaded {} contigs", mapper.contig_count());

    let r1_reader = fastq_reader(r1_path);
    let r2_reader = fastq_reader(r2_path);

    let mut batch = ReadBatch::with_capacity(batch_size);
    let mut total_reads = 0usize;
    let mut total_records = 0usize;
    let mut total_batches = 0usize;

    let start = Instant::now();

    let mut r1_iter = FastqIter::new(r1_reader);
    let mut r2_iter = FastqIter::new(r2_reader);

    loop {
        if total_reads >= max_reads {
            break;
        }

        // Fill batch
        batch.clear();
        for _ in 0..batch_size {
            if total_reads >= max_reads {
                break;
            }
            let Some(r1) = r1_iter.next() else { break };
            let Some(r2) = r2_iter.next() else { break };
            batch
                .push(
                    &r1.name,
                    ReadEnd { seq: &r1.seq, qual: &r1.qual },
                    ReadEnd { seq: &r2.seq, qual: &r2.qual },
                )
                .expect("failed to add read pair");
            total_reads += 1;
        }

        if batch.is_empty() {
            break;
        }

        let records = mapper.map_paired(&batch).unwrap_or_else(|e| {
            eprintln!("Mapping failed: {e}");
            std::process::exit(1);
        });
        total_records += records.len();
        total_batches += 1;

        if total_batches % 10 == 0 {
            eprintln!(
                "  batch {total_batches}: {total_reads} reads, {total_records} records, {:.1}s",
                start.elapsed().as_secs_f64()
            );
        }
    }

    let elapsed = start.elapsed().as_secs_f64();
    eprintln!("Done: {total_reads} reads in {total_batches} batches, {total_records} records");
    #[allow(clippy::cast_precision_loss)]
    let reads_per_sec = total_reads as f64 / elapsed;
    eprintln!("Time: {elapsed:.2}s ({reads_per_sec:.0} reads/sec)");
}

fn fastq_reader(path: &str) -> BufReader<Box<dyn std::io::Read>> {
    let file = File::open(path).unwrap_or_else(|e| {
        eprintln!("Cannot open {path}: {e}");
        std::process::exit(1);
    });
    let reader: Box<dyn std::io::Read> =
        if std::path::Path::new(path).extension().is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
        {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
    BufReader::with_capacity(256 * 1024, reader)
}

struct FastqRecord {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

struct FastqIter<R: BufRead> {
    reader: R,
    buf: String,
}

impl<R: BufRead> FastqIter<R> {
    fn new(reader: R) -> Self {
        Self { reader, buf: String::new() }
    }

    fn next(&mut self) -> Option<FastqRecord> {
        // Header line
        self.buf.clear();
        if self.reader.read_line(&mut self.buf).ok()? == 0 {
            return None;
        }
        let name = self.buf.trim_start_matches('@').split_whitespace().next()?.to_string();

        // Sequence
        self.buf.clear();
        if self.reader.read_line(&mut self.buf).ok()? == 0 {
            return None;
        }
        let seq = self.buf.trim_end().as_bytes().to_vec();

        // Plus line
        self.buf.clear();
        if self.reader.read_line(&mut self.buf).ok()? == 0 || !self.buf.starts_with('+') {
            return None;
        }

        // Quality
        self.buf.clear();
        if self.reader.read_line(&mut self.buf).ok()? == 0 {
            return None;
        }
        let qual = self.buf.trim_end().as_bytes().to_vec();

        assert_eq!(
            seq.len(),
            qual.len(),
            "FASTQ record '{name}': seq length {} != qual length {}",
            seq.len(),
            qual.len()
        );

        Some(FastqRecord { name, seq, qual })
    }
}
