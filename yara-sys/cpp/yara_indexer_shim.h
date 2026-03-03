// yara_indexer_shim.h — C-compatible API for the YARA FM index builder.
//
// This header declares an opaque handle to a YARA Indexer and a small set of
// functions for building an FM index from a FASTA file and querying the
// resulting contig metadata.

#ifndef YARA_INDEXER_SHIM_H_
#define YARA_INDEXER_SHIM_H_

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// Opaque handle
// ---------------------------------------------------------------------------

typedef struct YaraIndexerHandle YaraIndexerHandle;

// ---------------------------------------------------------------------------
// Configuration struct
// ---------------------------------------------------------------------------

typedef struct {
    const char* output_prefix;   // Index file prefix (e.g. "/tmp/myref")
    const char* tmp_dir;         // Temp directory for index construction (NULL = use output dir)
    int         verbose;         // 0, 1, 2
} YaraIndexerOptions;

// ---------------------------------------------------------------------------
// API functions
// ---------------------------------------------------------------------------

/// Build an FM index from a FASTA file.
/// Returns an opaque handle (for contig queries), or NULL on error.
/// error_buf is filled with a message on failure.
YaraIndexerHandle* yara_indexer_build(
    const char* fasta_path,
    const YaraIndexerOptions* opts,
    char* error_buf,
    size_t error_buf_len
);

/// Number of contigs in the built index.
size_t yara_indexer_contig_count(const YaraIndexerHandle* handle);

/// Name of the contig at `idx`.  Valid for the lifetime of the handle.
const char* yara_indexer_contig_name(const YaraIndexerHandle* handle, size_t idx);

/// Length of the contig at `idx`.
size_t yara_indexer_contig_length(const YaraIndexerHandle* handle, size_t idx);

/// Free the indexer handle and all associated memory.
void yara_indexer_close(YaraIndexerHandle* handle);

#ifdef __cplusplus
}
#endif

#endif // YARA_INDEXER_SHIM_H_
