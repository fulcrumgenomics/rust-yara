// yara_shim.cpp — C++ shim wrapping YARA's template-heavy mapper code.
//
// This file replicates the include chain from mapper.cpp, instantiates a
// subset of YARA's template configurations, and exposes C-callable functions
// that load an FM index, ingest reads from memory, run the mapping pipeline,
// and return alignment results as flat C structs instead of SAM text.

// ============================================================================
// Configuration
// ============================================================================

#define YARA_MAPPER
// YARA_LARGE_CONTIGS is defined via build.rs compiler flags.

// ============================================================================
// STL headers
// ============================================================================

#include <random>
#include <string>
#include <vector>
#include <memory>
#include <cstring>
#include <cmath>
#include <algorithm>

// ============================================================================
// SeqAn headers
// ============================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/parallel.h>

// ============================================================================
// YARA app headers (verbatim from seqan/apps/yara/)
// ============================================================================

// Forward declaration required by YARA headers before mapper.h defines it.
struct Options;

#include "basic_alphabet.h"
#include "file_pair.h"
#include "file_prefetched.h"
#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_tags.h"
#include "misc_types.h"
#include "index_fm.h"
#include "bits_reads.h"
#include "bits_hits.h"
#include "bits_context.h"
#include "bits_matches.h"
#include "bits_seeds.h"
#include "bits_bucket.h"
#include "find_verifier.h"
#include "find_extender.h"
#include "misc_options.h"
#include "mapper_collector.h"
#include "mapper_classifier.h"
#include "mapper_ranker.h"
#include "mapper_filter.h"
#include "mapper_extender.h"
#include "mapper_verifier.h"
#include "mapper_aligner.h"
#include "mapper_writer.h"
#include "mapper.h"

using namespace seqan2;

// ============================================================================
// Shim header
// ============================================================================

#include "yara_shim.h"

// ============================================================================
// Internal helpers
// ============================================================================

namespace {

// Convert a CIGAR string to a heap-allocated BAM-encoded uint32 array.
// BAM encoding: each op is (len << 4) | op_code
// Op codes: M=0, I=1, D=2, N=3, S=4, H=5, P=6, =7, X=8
static uint8_t cigarOpCode(char c) {
    switch (c) {
        case 'M': return 0;
        case 'I': return 1;
        case 'D': return 2;
        case 'N': return 3;
        case 'S': return 4;
        case 'H': return 5;
        case 'P': return 6;
        case '=': return 7;
        case 'X': return 8;
        default:  return 0;
    }
}

// ---------------------------------------------------------------------------
// FieldPools — bump allocator for seq/qual/cigar fields.
//
// Reduces per-record heap allocations to two bulk allocations per batch.
// XA strings remain individually heap-allocated (unpredictable sizes).
// The pool is owned by MapperInstance and reused across mapPaired calls.
// ---------------------------------------------------------------------------

struct FieldPools {
    char*     char_buf;
    size_t    char_cap;
    size_t    char_used;

    uint32_t* cigar_buf;
    size_t    cigar_cap;
    size_t    cigar_used;

    FieldPools() : char_buf(nullptr), char_cap(0), char_used(0),
                   cigar_buf(nullptr), cigar_cap(0), cigar_used(0) {}
    ~FieldPools() { delete[] char_buf; delete[] cigar_buf; }

    FieldPools(FieldPools const &) = delete;
    FieldPools& operator=(FieldPools const &) = delete;

    void prepare(size_t char_need, size_t cigar_need) {
        if (char_need > char_cap) {
            delete[] char_buf;
            char_buf = new char[char_need];
            char_cap = char_need;
        }
        char_used = 0;
        if (cigar_need > cigar_cap) {
            delete[] cigar_buf;
            cigar_buf = new uint32_t[cigar_need];
            cigar_cap = cigar_need;
        }
        cigar_used = 0;
    }

    char* alloc_chars(size_t n) {
        if (char_used + n > char_cap)
            throw std::runtime_error("FieldPools char buffer exhausted");
        char* p = char_buf + char_used;
        char_used += n;
        return p;
    }

    uint32_t* alloc_cigar(size_t n) {
        if (cigar_used + n > cigar_cap)
            throw std::runtime_error("FieldPools cigar buffer exhausted");
        uint32_t* p = cigar_buf + cigar_used;
        cigar_used += n;
        return p;
    }
};

template <typename TCigar>
static void encodeCigar(TCigar const & cigar, uint32_t** out_ops, uint32_t* out_len,
                        FieldPools& pools) {
    auto n = length(cigar);
    if (n == 0) {
        *out_ops = nullptr;
        *out_len = 0;
        return;
    }
    uint32_t* ops = pools.alloc_cigar(n);
    for (unsigned i = 0; i < n; ++i) {
        auto const & el = cigar[i];
        ops[i] = (static_cast<uint32_t>(el.count) << 4) | cigarOpCode(el.operation);
    }
    *out_ops = ops;
    *out_len = static_cast<uint32_t>(n);
}

// Copy a SeqAn Dna5Q string to a C string (ACGTN text).
// Pool overload: allocates from the provided FieldPools.
static const char BASES_TABLE[] = "ACGTN";

template <typename TSeq>
static char* seqToString(TSeq const & seq, FieldPools& pools) {
    auto n = length(seq);
    char* buf = pools.alloc_chars(n + 1);
    for (unsigned i = 0; i < n; ++i) {
        buf[i] = BASES_TABLE[std::min(static_cast<unsigned>(ordValue(seq[i])), 4u)];
    }
    buf[n] = '\0';
    return buf;
}

// Extract quality string from a Dna5Q sequence as phred+33.
// Pool overload: allocates from the provided FieldPools.
template <typename TSeq>
static char* qualToString(TSeq const & seq, FieldPools& pools) {
    auto n = length(seq);
    char* buf = pools.alloc_chars(n + 1);
    for (unsigned i = 0; i < n; ++i) {
        buf[i] = static_cast<char>(getQualityValue(seq[i]) + 33);
    }
    buf[n] = '\0';
    return buf;
}

// Duplicate a SeqAn CharString / std::string to a heap-allocated C string.
static char* dupString(char const* s, size_t len) {
    if (!s || len == 0) return nullptr;
    char* buf = new char[len + 1];
    std::memcpy(buf, s, len);
    buf[len] = '\0';
    return buf;
}

// Write an error message into the caller-provided buffer.
static void writeError(char* buf, size_t buf_len, const char* msg) {
    if (buf && buf_len > 0) {
        std::strncpy(buf, msg, buf_len - 1);
        buf[buf_len - 1] = '\0';
    }
}

// Free dynamically allocated fields from record structs.
// Pool-managed records (_pool_managed != 0): only xa is individually freed.
// Non-pool records (_pool_managed == 0): cigar, seq, qual, and xa are all freed.
static void freeRecordFields(YaraAlignmentRecord* records, size_t count) {
    if (!records) return;
    for (size_t i = 0; i < count; ++i) {
        if (!records[i]._pool_managed) {
            delete[] records[i].cigar;
            delete[] records[i].seq;
            delete[] records[i].qual;
        }
        delete[] records[i].xa;
        records[i].cigar = nullptr;
        records[i].seq = nullptr;
        records[i].qual = nullptr;
        records[i].xa = nullptr;
    }
}

// ============================================================================
// RecordCollector — replaces MatchesWriter, populating C structs instead of SAM.
// ============================================================================
//
// This is modeled on MatchesWriter in mapper_writer.h but writes to a flat
// array of YaraAlignmentRecord instead of a BAM file.

template <typename TSpec, typename Traits>
struct RecordCollector
{
    typedef typename Traits::TReads             TReads;
    typedef typename Traits::TMatches           TMatches;
    typedef typename Traits::TMatchesSet        TMatchesSet;
    typedef typename Traits::TMatchesViewView   TMatchesView;
    typedef typename Traits::TMatchesViewSet    TMatchesViewSet;
    typedef typename Traits::TMatchesProbs      TMatchesProbs;
    typedef typename Traits::TCigarsView        TCigarsView;
    typedef typename Traits::TCigarsSet         TCigarsSet;
    typedef typename Traits::TReadsContext      TReadsContext;

    // Output — out_count is a pointer so that copies made by iterate()
    // all increment the same counter (SeqAn's iterate passes functors by value).
    YaraAlignmentRecord*    out;
    size_t                  out_capacity;
    size_t*                 out_count_ptr;

    // Shared-memory read-only data.
    TMatchesViewSet const & matchesSet;
    TMatchesView const &    primaryMatches;
    TMatchesProbs const &   primaryMatchesProbs;
    TCigarsView const &     primaryCigars;
    TCigarsSet const &      cigarSet;
    TReadsContext const &   ctx;
    TReads const &          reads;
    Options const &         options;

    // Contig names for XA tag construction (TAG mode).
    std::vector<std::string> const * contigNames;

    // Pool allocator for seq/qual/cigar fields (owned by MapperInstance).
    FieldPools* pools;

    // Owned counter — lives on the original object, pointed to by copies.
    size_t                  out_count_storage;

    RecordCollector(YaraAlignmentRecord* out,
                    size_t out_capacity,
                    TMatchesViewSet const & matchesSet,
                    TMatchesView const & primaryMatches,
                    TMatchesProbs const & primaryMatchesProbs,
                    TCigarsView const & primaryCigars,
                    TCigarsSet const & cigarSet,
                    TReadsContext const & ctx,
                    TReads const & reads,
                    Options const & options,
                    std::vector<std::string> const * contigNames,
                    FieldPools* pools) :
        out(out),
        out_capacity(out_capacity),
        out_count_ptr(nullptr),
        matchesSet(matchesSet),
        primaryMatches(primaryMatches),
        primaryMatchesProbs(primaryMatchesProbs),
        primaryCigars(primaryCigars),
        cigarSet(cigarSet),
        ctx(ctx),
        reads(reads),
        options(options),
        contigNames(contigNames),
        pools(pools),
        out_count_storage(0)
    {
        out_count_ptr = &out_count_storage;
        // Process all primary matches (serially to avoid races on out_count).
        iterate(primaryMatches, *this, Standard(), Serial());
    }

    /// Number of records written.
    size_t out_count() const { return out_count_storage; }

    template <typename TIterator>
    void operator() (TIterator const & it)
    {
        _collectMatchesImpl(*this, it);
    }
};

// Forward declarations
template <typename TSpec, typename Traits, typename TMatchIt>
static void _collectMatchesImpl(RecordCollector<TSpec, Traits> & me, TMatchIt const & it);

template <typename TSpec, typename Traits, typename TReadId>
static void _collectUnmappedRead(RecordCollector<TSpec, Traits> & me, TReadId readId);

template <typename TSpec, typename Traits, typename TReadId, typename TMatch>
static void _collectMappedRead(RecordCollector<TSpec, Traits> & me, TReadId readId, TMatch const & primary);

// Allocate the next output slot, returning nullptr if full.
// The caller (mapPaired) zero-initializes the entire output buffer before use,
// so individual records do not need to be zeroed here.
template <typename TSpec, typename Traits>
static YaraAlignmentRecord* _nextRecord(RecordCollector<TSpec, Traits> & me) {
    if (*me.out_count_ptr >= me.out_capacity) return nullptr;
    return &me.out[(*me.out_count_ptr)++];
}

// Get the read pair index (0-based) from a readId.
// In paired-end mode: pairsCount = length(seqs)/4
// readId < pairsCount => first mate, pair_index = readId
// readId >= pairsCount => second mate, pair_index = readId - pairsCount
template <typename TReadSeqs, typename TReadId>
static uint32_t getPairIndex(TReadSeqs const & readSeqs, TReadId readId) {
    auto pairsCount = getPairsCount(readSeqs);
    if (static_cast<uint64_t>(readId) < pairsCount)
        return static_cast<uint32_t>(readId);
    return static_cast<uint32_t>(readId - pairsCount);
}

template <typename TReadSeqs, typename TReadId>
static uint8_t isRead1(TReadSeqs const & readSeqs, TReadId readId) {
    return static_cast<uint64_t>(readId) < getPairsCount(readSeqs) ? 1 : 0;
}

// Return BAM_FLAG_FIRST or BAM_FLAG_LAST based on mate position.
template <typename TReadSeqs, typename TReadId>
static uint16_t mateFlags(TReadSeqs const & readSeqs, TReadId readId) {
    return isFirstMate(readSeqs, readId) ? BAM_FLAG_FIRST : BAM_FLAG_LAST;
}

// Fill XA tag string for secondary matches.
// Format: contig_name,strand+pos,CIGAR,NM; (when alignSecondary is true)
// or:     contig_name,begin,end,strand,NM; (when alignSecondary is false)
template <typename TSpec, typename Traits, typename TMatches, typename TIter>
static char* _buildXa(RecordCollector<TSpec, Traits> & me, TMatches const & matches, TIter const & primaryIt) {
    typedef typename Value<TMatches const>::Type TMatch;
    std::string xa;
    xa.reserve(256);

    iterate(matches, [&](typename Iterator<TMatches const>::Type const & matchIt)
    {
        if (matchIt == primaryIt)
            return;
        TMatch const & match = *matchIt;

        // Contig identifier: use name if available, otherwise numeric ID
        auto contigId = getMember(match, ContigId());
        if (me.contigNames && contigId < me.contigNames->size()) {
            xa += (*me.contigNames)[contigId];
        } else {
            xa += std::to_string(contigId);
        }
        xa += ',';
        if (me.options.alignSecondary) {
            xa += onForwardStrand(match) ? '+' : '-';
            xa += std::to_string(getMember(match, ContigBegin()) + 1);
            xa += ',';
            // CIGAR for secondary
            auto const & cigar = me.cigarSet[getMember(match, ReadId())][position(matchIt, matches)];
            for (unsigned i = 0; i < length(cigar); ++i) {
                xa += std::to_string(cigar[i].count);
                xa += cigar[i].operation;
            }
            xa += ',';
            xa += std::to_string(getMember(match, Errors()));
            xa += ';';
        } else {
            xa += std::to_string(getMember(match, ContigBegin()) + 1);
            xa += ',';
            xa += std::to_string(getMember(match, ContigEnd()) + 1);
            xa += ',';
            xa += onForwardStrand(match) ? '+' : '-';
            xa += ',';
            xa += std::to_string(getMember(match, Errors()));
            xa += ';';
        }
    }, Standard(), Serial());

    if (xa.empty()) return nullptr;
    return dupString(xa.c_str(), xa.size());
}

// Dispatch: mapped vs unmapped
template <typename TSpec, typename Traits, typename TMatchIt>
static void _collectMatchesImpl(RecordCollector<TSpec, Traits> & me, TMatchIt const & it) {
    typedef typename Value<TMatchIt const>::Type TMatch;
    TMatch const & primary = value(it);
    if (isValid(primary))
        _collectMappedRead(me, position(it, me.primaryMatches), primary);
    else
        _collectUnmappedRead(me, position(it, me.primaryMatches));
}

// Write unmapped read pair entry (paired-end).
template <typename TSpec, typename Traits, typename TReadId>
static void _collectUnmappedRead(RecordCollector<TSpec, Traits> & me, TReadId readId) {
    auto* rec = _nextRecord(me);
    if (!rec) return;

    rec->read_pair_index = getPairIndex(me.reads.seqs, readId);
    rec->is_read1 = isRead1(me.reads.seqs, readId);
    rec->is_unmapped = 1;
    rec->flag = BAM_FLAG_UNMAPPED | BAM_FLAG_MULTIPLE | mateFlags(me.reads.seqs, readId);

    TReadId mateId = getMateId(me.reads.seqs, readId);
    if (!isMapped(me.ctx, mateId))
        rec->flag |= BAM_FLAG_NEXT_UNMAPPED;

    // Provide sequence and quality for unmapped reads.
    auto readSeqId = getFirstMateFwdSeqId(me.reads.seqs, readId);
    if (!isRead1(me.reads.seqs, readId))
        readSeqId = getSecondMateFwdSeqId(me.reads.seqs, getPairIndex(me.reads.seqs, readId));
    rec->seq = seqToString(me.reads.seqs[readSeqId], *me.pools);
    rec->qual = qualToString(me.reads.seqs[readSeqId], *me.pools);
    rec->seq_len = length(me.reads.seqs[readSeqId]);
    rec->_pool_managed = 1;
}

// Write mapped read (paired-end).
template <typename TSpec, typename Traits, typename TReadId, typename TMatch>
static void _collectMappedRead(RecordCollector<TSpec, Traits> & me, TReadId readId, TMatch const & primary) {
    typedef typename Traits::TMatches                           TMatches;
    typedef typename Size<TMatches>::Type                       TSize;
    typedef typename Iterator<TMatches const, Standard>::Type   TIter;

    auto* rec = _nextRecord(me);
    if (!rec) return;

    rec->read_pair_index = getPairIndex(me.reads.seqs, readId);
    rec->is_read1 = isRead1(me.reads.seqs, readId);
    rec->contig_id = getMember(primary, ContigId());
    rec->pos = getMember(primary, ContigBegin());
    rec->is_reverse = onReverseStrand(primary) ? 1 : 0;
    rec->is_secondary = 0;
    rec->is_unmapped = 0;
    rec->nm = getMember(primary, Errors());

    // SAM flags
    uint16_t flag = BAM_FLAG_MULTIPLE | mateFlags(me.reads.seqs, readId);
    if (onReverseStrand(primary)) flag |= BAM_FLAG_RC;

    // Mate info
    TReadId mateId = getMateId(me.reads.seqs, readId);
    if (!isMapped(me.ctx, mateId)) {
        flag |= BAM_FLAG_NEXT_UNMAPPED;
        // RNEXT/PNEXT = own position when mate unmapped
        rec->mate_contig_id = rec->contig_id;
        rec->mate_pos = rec->pos;
    } else {
        TMatch const & mate = me.primaryMatches[mateId];
        rec->mate_contig_id = getMember(mate, ContigId());
        rec->mate_pos = getMember(mate, ContigBegin());
        if (onReverseStrand(mate)) flag |= BAM_FLAG_NEXT_RC;

        // Proper pair and template length
        if (isPaired(me.ctx, readId)) {
            if (orientationProper(primary, mate))
                flag |= BAM_FLAG_ALL_PROPER;

            if (getMember(primary, ContigId()) == getMember(mate, ContigId())) {
                if (getMember(primary, ContigBegin()) < getMember(mate, ContigBegin()))
                    rec->tlen = getMember(mate, ContigEnd()) - getMember(primary, ContigBegin());
                else
                    rec->tlen = getMember(mate, ContigBegin()) - getMember(primary, ContigEnd());
            }
        }
    }
    rec->flag = flag;

    // Match counts and MAPQ
    TMatches const & matches = me.matchesSet[readId];
    TSize bestCount = countMatchesInBestStratum(matches);
    TSize subCount = length(matches) - bestCount;
    rec->x0 = static_cast<uint16_t>(std::min(bestCount, static_cast<TSize>(65535)));
    rec->x1 = static_cast<uint16_t>(std::min(subCount, static_cast<TSize>(65535)));

    if (isPaired(me.ctx, readId)) {
        rec->mapq = getMapq(me.primaryMatchesProbs[readId]);
    } else {
        double errorRate = getErrorRate(primary, me.reads.seqs);
        double prob = getMatchProb(errorRate, errorRate, bestCount, subCount);
        rec->mapq = getMapq(prob);
    }

    // CIGAR
    auto const & cigar = me.primaryCigars[getMember(primary, ReadId())];
    encodeCigar(cigar, &rec->cigar, &rec->cigar_len, *me.pools);

    // Sequence and quality for primary
    auto readSeqId = getReadSeqId(primary, me.reads.seqs);
    rec->seq = seqToString(me.reads.seqs[readSeqId], *me.pools);
    rec->qual = qualToString(me.reads.seqs[readSeqId], *me.pools);
    rec->seq_len = length(me.reads.seqs[readSeqId]);
    rec->_pool_managed = 1;

    // XA tag for secondary matches
    TIter primaryIt = findMatch(matches, primary);
    if (me.options.secondaryMatches == TAG) {
        rec->xa = _buildXa(me, matches, primaryIt);
    }

    // RECORD mode: emit separate secondary records
    if (me.options.secondaryMatches == RECORD) {
        iterate(matches, [&](typename Iterator<TMatches const>::Type matchIt)
        {
            if (matchIt == primaryIt) return;
            TMatch const & match = *matchIt;

            auto* srec = _nextRecord(me);
            if (!srec) return;

            srec->read_pair_index = rec->read_pair_index;
            srec->is_read1 = rec->is_read1;
            srec->contig_id = getMember(match, ContigId());
            srec->pos = getMember(match, ContigBegin());
            srec->is_reverse = onReverseStrand(match) ? 1 : 0;
            srec->is_secondary = 1;
            srec->nm = getMember(match, Errors());

            uint16_t sflag = BAM_FLAG_MULTIPLE | BAM_FLAG_SECONDARY | mateFlags(me.reads.seqs, readId);
            if (onReverseStrand(match)) sflag |= BAM_FLAG_RC;
            srec->flag = sflag;

            // CIGAR for secondary (if align_secondary)
            if (me.options.alignSecondary) {
                auto const & scigar = me.cigarSet[getMember(match, ReadId())][position(matchIt, matches)];
                encodeCigar(scigar, &srec->cigar, &srec->cigar_len, *me.pools);
            }

            // Secondaries don't get seq/qual (SEQ=*, QUAL=*)
            srec->seq = nullptr;
            srec->qual = nullptr;
            srec->seq_len = 0;
            srec->_pool_managed = 1;
        }, Standard(), Serial());
    }
}

// ============================================================================
// Read ingestion helper
// ============================================================================

// Append reads from C string arrays into a SeqAn StringSet (Dna5Q).
template <typename TReadSeq, typename TReadSeqs>
static void ingestReads(
    TReadSeqs & seqs,
    size_t count,
    const char* const* sequences,
    const char* const* qualities
) {
    for (size_t i = 0; i < count; ++i) {
        auto len = std::strlen(sequences[i]);
        TReadSeq seq;
        resize(seq, len);
        for (size_t j = 0; j < len; ++j) {
            switch (sequences[i][j]) {
                case 'A': case 'a': seq[j] = Dna5Q('A'); break;
                case 'C': case 'c': seq[j] = Dna5Q('C'); break;
                case 'G': case 'g': seq[j] = Dna5Q('G'); break;
                case 'T': case 't': seq[j] = Dna5Q('T'); break;
                default:            seq[j] = Dna5Q('N'); break;
            }
            assignQualityValue(seq[j], qualities[i][j] - 33);
        }
        appendValue(seqs, seq);
    }
}

// ============================================================================
// Type-erased mapper handle
// ============================================================================

// Virtual base class so the handle can hold any template instantiation.

struct MapperBase {
    virtual ~MapperBase() = default;
    virtual int64_t mapPaired(
        YaraReadBatch const* reads,
        YaraAlignmentRecord* out,
        size_t out_capacity,
        char* error_buf,
        size_t error_buf_len
    ) = 0;
    virtual size_t contigCount() const = 0;
    virtual const char* contigName(size_t idx) const = 0;
    virtual size_t contigLength(size_t idx) const = 0;
};

template <typename TConfig>
struct MapperInstance : MapperBase {
    Options options;
    Mapper<void, TConfig> mapper;

    // Cached contig name C-strings for the lifetime of the handle.
    std::vector<std::string> contigNameStrings;

    // Pool allocator for seq/qual/cigar fields, reused across mapPaired calls.
    FieldPools pools;

    MapperInstance(Options const & opts) : options(opts), mapper(options) {}

    void loadIndex() {
        configureThreads(mapper);
        loadContigs(mapper);
        loadContigsIndex(mapper);

        // Cache contig names as std::string for C API access.
        // contigs.names is a ConcatDirect StringSet, so elements are Segments.
        contigNameStrings.resize(length(mapper.contigs.names));
        for (unsigned i = 0; i < length(mapper.contigs.names); ++i) {
            auto const & name = mapper.contigs.names[i];
            std::string s;
            s.reserve(length(name));
            for (unsigned j = 0; j < length(name); ++j)
                s += static_cast<char>(name[j]);
            contigNameStrings[i] = std::move(s);
        }
    }

    size_t contigCount() const override {
        return length(mapper.contigs.seqs);
    }

    const char* contigName(size_t idx) const override {
        if (idx >= contigNameStrings.size()) return "";
        return contigNameStrings[idx].c_str();
    }

    size_t contigLength(size_t idx) const override {
        if (idx >= length(mapper.contigs.seqs)) return 0;
        return length(mapper.contigs.seqs[idx]);
    }

    int64_t mapPaired(
        YaraReadBatch const* batch,
        YaraAlignmentRecord* out,
        size_t out_capacity,
        char* error_buf,
        size_t error_buf_len
    ) override {
        typedef MapperTraits<void, TConfig> TTraits;
        typedef typename TTraits::TMatch       TMatch;
        typedef typename TTraits::TReads::TSeq TReadSeq;

        // Zero-initialize the output buffer so all pointer fields start as
        // null.  This is required for safe cleanup in the error path.
        std::memset(out, 0, out_capacity * sizeof(YaraAlignmentRecord));

        try {
            // 1. Populate reads store from the C batch.
            clear(mapper.reads);
            auto & seqs = mapper.reads.seqs;
            auto & names = mapper.reads.names;

            size_t n = batch->count;
            if (n == 0) return 0;

            reserve(seqs, 4 * n, Exact());
            reserve(concat(seqs), 4 * n * 200, Exact());
            reserve(names, n, Exact());

            // Append R1 forward sequences (indices 0..n-1)
            ingestReads<TReadSeq>(seqs, n, batch->r1_seqs, batch->r1_quals);

            // Append R2 forward sequences (indices n..2n-1)
            ingestReads<TReadSeq>(seqs, n, batch->r2_seqs, batch->r2_quals);

            // Names
            for (size_t i = 0; i < n; ++i) {
                CharString name(batch->names[i]);
                appendValue(names, name);
            }

            // Check max read length
            if (maxLength(seqs, typename TConfig::TThreading()) > MemberLimits<TMatch, ReadSize>::VALUE) {
                writeError(error_buf, error_buf_len, "Maximum read length exceeded.");
                return -1;
            }

            // Append reverse complements (layout: R1fwd, R2fwd, R1rc, R2rc)
            appendReverseComplement(mapper.reads);

            // 2. Run mapping pipeline + collect results.
            return mapReadsAndCollect(mapper, seqs, out, out_capacity);

        } catch (std::exception const & e) {
            // Free any heap-allocated fields in partially written records
            // before returning an error, so the Rust caller doesn't need to.
            freeRecordFields(out, out_capacity);
            writeError(error_buf, error_buf_len, e.what());
            return -1;
        }
    }

private:
    // Replicate _mapReadsImpl but collect to C structs instead of writing SAM.
    template <typename TReadSeqs>
    int64_t mapReadsAndCollect(
        Mapper<void, TConfig> & me,
        TReadSeqs & readSeqs,
        YaraAlignmentRecord* out,
        size_t out_capacity
    ) {
        typedef MapperTraits<void, TConfig> TTraits;

        initReadsContext(me, readSeqs);
        initSeeds(me, readSeqs);

        collectSeeds<0>(me, readSeqs);
        findSeeds<0>(me, 0);
        classifyReads(me);
        collectSeeds<1>(me, readSeqs);
        collectSeeds<2>(me, readSeqs);
        findSeeds<0>(me, 1);
        findSeeds<0>(me, 2);
        rankSeeds(me);
        reserveMatches(me);
        extendHits<0>(me, 0);
        extendHits<0>(me, 1);
        extendHits<0>(me, 2);
        clearSeeds(me);
        clearHits(me);

        initSeeds(me, readSeqs);
        collectSeeds<1>(me, readSeqs);
        findSeeds<1>(me, 1);
        collectSeeds<2>(me, readSeqs);
        findSeeds<1>(me, 2);
        rankSeeds(me);
        extendHits<1>(me, 1);
        extendHits<1>(me, 2);
        clearSeeds(me);
        clearHits(me);

        if (me.options.sensitivity > LOW) {
            initSeeds(me, readSeqs);
            collectSeeds<2>(me, readSeqs);
            findSeeds<2>(me, 2);
            rankSeeds(me);
            extendHits<2>(me, 2);
            clearHits(me);
            clearSeeds(me);
        }

        aggregateMatches(me, readSeqs);
        rankMatches(me, readSeqs);
        if (me.options.verifyMatches)
            verifyMatches(me);
        alignMatches(me);

        // Prepare pools for seq/qual/cigar allocation.
        // Exact sizing for seq+qual (forward reads only); generous for cigar.
        {
            size_t n_fwd = length(readSeqs) / 2;
            size_t char_need = 0;
            for (size_t i = 0; i < n_fwd; ++i) {
                char_need += (length(readSeqs[i]) + 1) * 2; // seq + qual per read
            }
            size_t cigar_need = out_capacity * 64; // 64 ops per record (generous)
            pools.prepare(char_need, cigar_need);
        }

        // Collect results into C structs instead of writing SAM.
        RecordCollector<void, TTraits> collector(
            out, out_capacity,
            me.matchesSet,
            me.primaryMatches, me.primaryMatchesProbs,
            me.primaryCigars, me.cigarsSet,
            me.ctx, me.reads,
            me.options,
            &contigNameStrings,
            &pools
        );

        int64_t count = static_cast<int64_t>(collector.out_count());

        clearMatches(me);
        clearAlignments(me);

        return count;
    }
};

} // anonymous namespace

// ============================================================================
// The opaque handle wraps a MapperBase pointer.
// ============================================================================

struct YaraMapperHandle {
    std::unique_ptr<MapperBase> impl;
};

// ============================================================================
// Template instantiation factory
// ============================================================================

// We instantiate a subset of template configs needed for the fghla use case
// and a general-purpose set.  The contig-size types are selected at runtime
// based on the .txt.size file, just like mapper.cpp does.

namespace {

template <typename TContigsSize, typename TContigsLen, typename TContigsSum,
          typename TSeedsDistance>
static MapperBase* createMapper(Options const & options) {
    typedef ReadMapperConfig<Parallel, PairedEnd, TSeedsDistance,
                             TContigsSize, TContigsLen, TContigsSum> TConfig;
    auto* inst = new MapperInstance<TConfig>(options);
    inst->loadIndex();
    return inst;
}

template <typename TContigsSize, typename TContigsLen, typename TContigsSum>
static MapperBase* createMapperBySensitivity(Options const & options) {
    if (options.sensitivity == FULL)
        return createMapper<TContigsSize, TContigsLen, TContigsSum, EditDistance>(options);
    else
        return createMapper<TContigsSize, TContigsLen, TContigsSum, HammingDistance>(options);
}

template <typename TContigsSize, typename TContigsLen>
static MapperBase* createMapperByContigsSum(Options const & options) {
    if (options.contigsSum <= std::numeric_limits<uint32_t>::max())
        return createMapperBySensitivity<TContigsSize, TContigsLen, uint32_t>(options);
    else
        return createMapperBySensitivity<TContigsSize, TContigsLen, uint64_t>(options);
}

template <typename TContigsSize>
static MapperBase* createMapperByContigsLen(Options const & options) {
    if (options.contigsMaxLength <= std::numeric_limits<uint32_t>::max())
        return createMapperByContigsSum<TContigsSize, uint32_t>(options);
    else
        return createMapperByContigsSum<TContigsSize, uint64_t>(options);
}

static MapperBase* createMapperByContigsSize(Options const & options) {
    if (options.contigsSize <= std::numeric_limits<uint8_t>::max())
        return createMapperByContigsLen<uint8_t>(options);
    else if (options.contigsSize <= std::numeric_limits<uint16_t>::max())
        return createMapperByContigsLen<uint16_t>(options);
    else
        return createMapperByContigsLen<uint32_t>(options);
}

} // anonymous namespace

// ============================================================================
// C API implementation
// ============================================================================

extern "C" {

YaraMapperHandle* yara_mapper_open(
    const char* index_prefix,
    const YaraMapperOptions* opts,
    char* error_buf,
    size_t error_buf_len
) {
    try {
        // Translate C options to YARA Options struct.
        Options options;
        options.contigsIndexFile = index_prefix;
        options.errorRate = opts->error_rate;
        options.strataRate = opts->strata_rate;
        options.strataCount = (opts->strata_count < 0) ? static_cast<uint64_t>(-1) : opts->strata_count;
        options.sensitivity = static_cast<Sensitivity>(opts->sensitivity);
        options.threadsCount = std::max(opts->threads, 1);
        options.secondaryMatches = static_cast<SecondaryAlignments>(opts->secondary_mode);
        options.alignSecondary = opts->align_secondary != 0;
        options.verifyMatches = opts->verify_matches != 0;
        options.verbose = opts->verbose;
        options.singleEnd = false; // always paired-end
        options.libraryLength = opts->library_length;
        options.libraryDev = opts->library_dev;

        // Load contig limits from the .txt.size file.
        if (!openContigsLimits(options)) {
            writeError(error_buf, error_buf_len, "Error while opening reference index size file.");
            return nullptr;
        }

        // Create the mapper with the right template instantiation.
        auto* impl = createMapperByContigsSize(options);

        auto* handle = new YaraMapperHandle();
        handle->impl.reset(impl);
        return handle;

    } catch (std::exception const & e) {
        writeError(error_buf, error_buf_len, e.what());
        return nullptr;
    }
}

int64_t yara_mapper_map_paired(
    YaraMapperHandle* handle,
    const YaraReadBatch* reads,
    YaraAlignmentRecord* out,
    size_t out_capacity,
    char* error_buf,
    size_t error_buf_len
) {
    if (!handle || !handle->impl) {
        writeError(error_buf, error_buf_len, "Invalid mapper handle.");
        return -1;
    }
    return handle->impl->mapPaired(reads, out, out_capacity, error_buf, error_buf_len);
}

size_t yara_mapper_contig_count(const YaraMapperHandle* handle) {
    if (!handle || !handle->impl) return 0;
    return handle->impl->contigCount();
}

const char* yara_mapper_contig_name(const YaraMapperHandle* handle, size_t idx) {
    if (!handle || !handle->impl) return "";
    return handle->impl->contigName(idx);
}

size_t yara_mapper_contig_length(const YaraMapperHandle* handle, size_t idx) {
    if (!handle || !handle->impl) return 0;
    return handle->impl->contigLength(idx);
}

void yara_mapper_close(YaraMapperHandle* handle) {
    delete handle;
}

void yara_mapper_free_records(YaraAlignmentRecord* records, size_t count) {
    freeRecordFields(records, count);
}

void yara_mapper_free_record(YaraAlignmentRecord* record) {
    if (record) freeRecordFields(record, 1);
}

} // extern "C"
