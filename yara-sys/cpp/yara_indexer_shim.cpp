// yara_indexer_shim.cpp — C++ shim wrapping YARA's template-heavy indexer code.
//
// This file replicates the include chain from indexer.cpp, instantiates the
// template configurations needed for FM index building, and exposes C-callable
// functions that build an index from a FASTA file and return contig metadata.

// ============================================================================
// Configuration
// ============================================================================

#define YARA_INDEXER
// YARA_LARGE_CONTIGS is defined via build.rs compiler flags.

// ============================================================================
// STL headers
// ============================================================================

#include <string>
#include <vector>
#include <memory>
#include <cstring>
#include <cstdlib>

// ============================================================================
// SeqAn headers
// ============================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

// ============================================================================
// YARA app headers (verbatim from seqan/apps/yara/)
// ============================================================================

// Forward declaration required by YARA headers before we define our Options.
struct Options;

#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_tags.h"
#include "misc_types.h"
#include "bits_matches.h"
#include "misc_options.h"
#include "index_fm.h"

using namespace seqan2;

// ============================================================================
// Shim header
// ============================================================================

#include "yara_indexer_shim.h"

// ============================================================================
// Internal Options struct
// ============================================================================

// Minimal Options struct compatible with the TOptions templates in
// misc_options.h.  Needs contigsIndexFile, contigsSize, contigsMaxLength,
// contigsSum for saveContigsLimits/setContigsLimits/factory dispatch.
struct Options {
    CharString contigsIndexFile;
    uint64_t contigsSize;
    uint64_t contigsMaxLength;
    uint64_t contigsSum;
    bool verbose;

    Options() :
        contigsSize(),
        contigsMaxLength(),
        contigsSum(),
        verbose(false)
    {}
};

// ============================================================================
// Internal helpers
// ============================================================================

namespace {

// Write an error message into the caller-provided buffer.
static void writeError(char* buf, size_t buf_len, const char* msg) {
    if (buf && buf_len > 0) {
        std::strncpy(buf, msg, buf_len - 1);
        buf[buf_len - 1] = '\0';
    }
}

// RAII guard for saving and restoring the TMPDIR environment variable.
// SeqAn uses TMPDIR for intermediate files during FM index construction.
// Since setenv is a global side-effect, we save the original value on
// construction and restore it on destruction (including exception unwind).
struct TmpdirGuard {
    std::string saved;
    bool had_value;

    TmpdirGuard() : had_value(false) {
        const char* old = std::getenv("TMPDIR");
        if (old) {
            saved = old;
            had_value = true;
        }
    }

    ~TmpdirGuard() {
        if (had_value) {
            CharString restore(saved);
            setEnv("TMPDIR", restore);
        }
    }

    TmpdirGuard(TmpdirGuard const &) = delete;
    TmpdirGuard& operator=(TmpdirGuard const &) = delete;
};

// ============================================================================
// Type-erased indexer handle
// ============================================================================

struct IndexerBase {
    virtual ~IndexerBase() = default;
    virtual size_t contigCount() const = 0;
    virtual const char* contigName(size_t idx) const = 0;
    virtual size_t contigLength(size_t idx) const = 0;
};

template <typename TContigsSize, typename TContigsLen, typename TContigsSum>
struct IndexerInstance : IndexerBase {
    // Cached contig metadata for queries after building.
    std::vector<std::string> contigNames;
    std::vector<size_t> contigLengths;

    IndexerInstance() = default;

    void build(SeqStore<void, YaraContigsConfig<>> & contigs, Options const & options) {
        typedef YaraFMConfig<TContigsSize, TContigsLen, TContigsSum>    TIndexConfig;
        typedef FMIndex<void, TIndexConfig>                             TIndexSpec;
        typedef Index<typename TIndexConfig::Text, TIndexSpec>          TIndex;

        // Cache contig metadata before index building clears them.
        auto n = length(contigs.seqs);
        contigNames.resize(n);
        contigLengths.resize(n);
        for (unsigned i = 0; i < n; ++i) {
            auto const & name = contigs.names[i];
            std::string s;
            s.reserve(length(name));
            for (unsigned j = 0; j < length(name); ++j)
                s += static_cast<char>(name[j]);
            contigNames[i] = std::move(s);
            contigLengths[i] = length(contigs.seqs[i]);
        }

        // Randomly replace Ns with A/C/G/T (seed 0xDEADBEEF).
        randomizeNs(contigs);

        // IndexFM is built on the reversed contigs.
        reverse(contigs);

        TIndex index;

        // Copy contigs into index text (Dna5 -> Dna conversion).
        setValue(index.text, contigs.seqs);

        // Clear contigs to free memory.
        clear(contigs);
        shrinkToFit(contigs);

        try {
            // Iterator instantiation triggers FM index construction.
            typename Iterator<TIndex, TopDown<>>::Type it(index);
            ignoreUnusedVariableWarning(it);
        }
        catch (BadAlloc const & /* e */) {
            throw RuntimeError("Insufficient memory to index the reference.");
        }
        catch (IOError const & /* e */) {
            throw RuntimeError("Insufficient disk space to index the reference. "
                               "Specify a bigger temporary folder using tmp_dir.");
        }

        // Save FM index files (.sa, .lf).
        if (!save(index, toCString(options.contigsIndexFile)))
            throw RuntimeError("Error while saving the reference index file.");
    }

    size_t contigCount() const override {
        return contigNames.size();
    }

    const char* contigName(size_t idx) const override {
        if (idx >= contigNames.size()) return "";
        return contigNames[idx].c_str();
    }

    size_t contigLength(size_t idx) const override {
        if (idx >= contigLengths.size()) return 0;
        return contigLengths[idx];
    }
};

// ============================================================================
// Template instantiation factory
// ============================================================================

template <typename TContigsSize, typename TContigsLen, typename TContigsSum>
static IndexerBase* createIndexer(
    SeqStore<void, YaraContigsConfig<>> & contigs,
    Options const & options
) {
    auto* inst = new IndexerInstance<TContigsSize, TContigsLen, TContigsSum>();
    inst->build(contigs, options);
    return inst;
}

template <typename TContigsSize, typename TContigsLen>
static IndexerBase* createIndexerByContigsSum(
    SeqStore<void, YaraContigsConfig<>> & contigs,
    Options const & options
) {
    if (options.contigsSum <= std::numeric_limits<uint32_t>::max())
        return createIndexer<TContigsSize, TContigsLen, uint32_t>(contigs, options);
    else
        return createIndexer<TContigsSize, TContigsLen, uint64_t>(contigs, options);
}

template <typename TContigsSize>
static IndexerBase* createIndexerByContigsLen(
    SeqStore<void, YaraContigsConfig<>> & contigs,
    Options const & options
) {
    if (options.contigsMaxLength <= std::numeric_limits<uint32_t>::max())
        return createIndexerByContigsSum<TContigsSize, uint32_t>(contigs, options);
    else
        return createIndexerByContigsSum<TContigsSize, uint64_t>(contigs, options);
}

static IndexerBase* createIndexerByContigsSize(
    SeqStore<void, YaraContigsConfig<>> & contigs,
    Options const & options
) {
    if (options.contigsSize <= std::numeric_limits<uint8_t>::max())
        return createIndexerByContigsLen<uint8_t>(contigs, options);
    else if (options.contigsSize <= std::numeric_limits<uint16_t>::max())
        return createIndexerByContigsLen<uint16_t>(contigs, options);
    else
        return createIndexerByContigsLen<uint32_t>(contigs, options);
}

} // anonymous namespace

// ============================================================================
// The opaque handle wraps an IndexerBase pointer.
// ============================================================================

struct YaraIndexerHandle {
    std::unique_ptr<IndexerBase> impl;
};

// ============================================================================
// C API implementation
// ============================================================================

extern "C" {

YaraIndexerHandle* yara_indexer_build(
    const char* fasta_path,
    const YaraIndexerOptions* opts,
    char* error_buf,
    size_t error_buf_len
) {
    try {
        Options options;
        options.contigsIndexFile = opts->output_prefix;
        options.verbose = opts->verbose > 0;

        // Save and restore TMPDIR via RAII (restored on success, exception,
        // or early return).
        TmpdirGuard tmpdirGuard;

        if (opts->tmp_dir != nullptr) {
            CharString tmpDir(opts->tmp_dir);
            setEnv("TMPDIR", tmpDir);
        } else {
            // Default: use the output directory as tmpdir.
            CharString indexFile(opts->output_prefix);
            CharString tmpDir = getPath(indexFile);
            if (empty(tmpDir))
                getCwd(tmpDir);
            setEnv("TMPDIR", tmpDir);
        }

        // 1. Load FASTA into SeqStore.
        typedef SeqStore<void, YaraContigsConfig<>> TContigs;
        TContigs contigs;
        SeqFileIn contigsFile;

        if (!open(contigsFile, fasta_path))
            throw RuntimeError("Error while opening the reference file.");

        try {
            readRecords(contigs, contigsFile);
            trimSeqNames(contigs);
        }
        catch (BadAlloc const & /* e */) {
            throw RuntimeError("Insufficient memory to load the reference.");
        }

        // 2. Compute contig limits.
        setContigsLimits(options, contigs.seqs);

        // 3. Save contigs (.txt, .rid) and limits (.txt.size).
        if (!saveContigsLimits(options) || !save(contigs, toCString(options.contigsIndexFile)))
            throw RuntimeError("Error while saving the reference.");

        // 4. Dispatch to factory chain: build and save FM index.
        auto* impl = createIndexerByContigsSize(contigs, options);

        auto* handle = new YaraIndexerHandle();
        handle->impl.reset(impl);
        return handle;

    } catch (std::exception const & e) {
        writeError(error_buf, error_buf_len, e.what());
        return nullptr;
    } catch (...) {
        writeError(error_buf, error_buf_len, "Unknown C++ exception during index build.");
        return nullptr;
    }
}

size_t yara_indexer_contig_count(const YaraIndexerHandle* handle) {
    if (!handle || !handle->impl) return 0;
    return handle->impl->contigCount();
}

const char* yara_indexer_contig_name(const YaraIndexerHandle* handle, size_t idx) {
    if (!handle || !handle->impl) return "";
    return handle->impl->contigName(idx);
}

size_t yara_indexer_contig_length(const YaraIndexerHandle* handle, size_t idx) {
    if (!handle || !handle->impl) return 0;
    return handle->impl->contigLength(idx);
}

void yara_indexer_close(YaraIndexerHandle* handle) {
    delete handle;
}

} // extern "C"
