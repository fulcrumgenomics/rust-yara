// yara_shim.h — C-compatible API for the YARA read mapper.
//
// This header declares an opaque handle to a YARA Mapper and a small set of
// functions for loading an FM index, mapping paired-end reads from memory, and
// retrieving alignment results as flat C structs (no SAM serialization).

#ifndef YARA_SHIM_H_
#define YARA_SHIM_H_

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// Opaque handle
// ---------------------------------------------------------------------------

typedef struct YaraMapperHandle YaraMapperHandle;

// ---------------------------------------------------------------------------
// Configuration structs
// ---------------------------------------------------------------------------

typedef struct {
    float error_rate;      // -e (0.05 = 5%)
    float strata_rate;     // -s (0.00 = 0%)
    int   strata_count;    // -sc (alternative to strata_rate), -1 = unused
    int   sensitivity;     // 0=LOW, 1=HIGH, 2=FULL
    int   threads;         // -t
    int   secondary_mode;  // 0=TAG, 1=RECORD, 2=OMIT
    int   align_secondary; // -as flag
    int   verify_matches;  // !-ni flag
    int   verbose;         // 0, 1, 2
    int   library_length;  // -ll (0 = auto-estimate)
    int   library_dev;     // -ld (0 = auto-estimate)
} YaraMapperOptions;

typedef struct {
    const char* const* names;       // read names array
    const char* const* r1_seqs;     // R1 sequences (ACGTN text)
    const char* const* r1_quals;    // R1 phred+33 quality strings
    const char* const* r2_seqs;     // R2 sequences
    const char* const* r2_quals;    // R2 phred+33 quality strings
    size_t count;                   // number of read pairs
} YaraReadBatch;

typedef struct {
    // Read identification
    uint32_t read_pair_index;   // index into input batch (0-based)
    uint8_t  is_read1;          // 1=R1, 0=R2

    // Alignment position
    uint32_t contig_id;         // reference contig index
    uint32_t pos;               // 0-based leftmost position
    uint8_t  is_reverse;        // mapped to reverse strand
    uint8_t  is_secondary;      // secondary alignment flag
    uint8_t  is_unmapped;       // unmapped flag

    // Alignment quality
    uint8_t  mapq;              // mapping quality
    uint8_t  nm;                // edit distance (NM tag)
    uint16_t x0;                // co-optimal alignment count (X0)
    uint16_t x1;                // sub-optimal alignment count (X1)

    // Mate info
    uint32_t mate_contig_id;
    uint32_t mate_pos;
    int32_t  tlen;              // template length
    uint16_t flag;              // full SAM flag

    // CIGAR (BAM-encoded uint32 array: op | len<<4)
    uint32_t* cigar;            // heap-allocated; freed by yara_mapper_free_records
    uint32_t cigar_len;         // number of CIGAR operations

    // Sequence and quality (NULL for secondaries / unmapped without seq)
    const char* seq;            // heap-allocated; freed by yara_mapper_free_records
    const char* qual;           // heap-allocated; freed by yara_mapper_free_records
    uint32_t seq_len;

    // XA tag (when secondary_mode=TAG, populated for primary records)
    const char* xa;             // heap-allocated; freed by yara_mapper_free_records
} YaraAlignmentRecord;

// ---------------------------------------------------------------------------
// API functions
// ---------------------------------------------------------------------------

/// Create mapper, load index.  Returns NULL on error.
/// error_buf is filled with a message on failure.
YaraMapperHandle* yara_mapper_open(
    const char* index_prefix,
    const YaraMapperOptions* opts,
    char* error_buf,
    size_t error_buf_len
);

/// Map a batch of paired-end reads.
/// Returns number of output records written to `out`, or -1 on error.
int64_t yara_mapper_map_paired(
    YaraMapperHandle* handle,
    const YaraReadBatch* reads,
    YaraAlignmentRecord* out,
    size_t out_capacity,
    char* error_buf,
    size_t error_buf_len
);

/// Number of reference contigs in the loaded index.
size_t yara_mapper_contig_count(const YaraMapperHandle* handle);

/// Name of the contig at `idx`.  Valid for the lifetime of the handle.
const char* yara_mapper_contig_name(const YaraMapperHandle* handle, size_t idx);

/// Length of the contig at `idx`.
size_t yara_mapper_contig_length(const YaraMapperHandle* handle, size_t idx);

/// Free the mapper and all associated memory.
void yara_mapper_close(YaraMapperHandle* handle);

/// Free heap-allocated fields in a batch of records (cigar, seq, qual, xa).
void yara_mapper_free_records(YaraAlignmentRecord* records, size_t count);

/// Free heap-allocated fields in a single record (cigar, seq, qual, xa).
void yara_mapper_free_record(YaraAlignmentRecord* record);

#ifdef __cplusplus
}
#endif

#endif // YARA_SHIM_H_
