# Yara Performance Optimization Implementation Plan

This document outlines a comprehensive plan for optimizing Yara's performance while maintaining functional output compatibility. Optimizations are grouped by priority and dependency.

## Overview

**Goal:** Improve throughput for HLA calling workflows (e.g., OptiType) without changing output format or alignment results.

**Estimated Total Effort:** 4-6 weeks for full implementation

**Testing Strategy:** All changes must pass existing tests in `tests/run_tests.py` with byte-identical SAM output (excluding @PG headers).

---

## Phase 1: Quick Wins (Week 1)

### 1.1 Increase Output Buffer Size

**File:** `mapper_writer.h:657`

**Current:**
```cpp
if (length(me.recordBuffer) > Power<2, 16>::VALUE)  // 64KB
```

**Change to:**
```cpp
if (length(me.recordBuffer) > Power<2, 20>::VALUE)  // 1MB
```

**Rationale:** Reduces critical section contention by 16x fewer flushes.

**Risk:** Low - only affects buffering, not output content.

**Validation:**
- Run with `-v 2` to compare I/O timing
- Verify identical output with `diff`

---

### 1.2 Pre-compute Sort Key Shift Constants

**File:** `bits_matches.h:505-528`

**Current:** Shift amounts computed inline in hot path
```cpp
return ((uint64_t)getMember(me, ContigId()) << (1 + MemberBits<TMatch, ContigSize>::VALUE + MemberBits<TMatch, Errors>::VALUE)) | ...
```

**New:** Add compile-time constant helper struct before `getSortKey` functions (~line 500):

```cpp
template <typename TSpec>
struct SortKeyShifts
{
    static constexpr unsigned ERRORS_BITS = MemberBits<Match<TSpec>, Errors>::VALUE;
    static constexpr unsigned CONTIG_SIZE_BITS = MemberBits<Match<TSpec>, ContigSize>::VALUE;

    static constexpr unsigned ERRORS_SHIFT = 0;
    static constexpr unsigned CONTIG_POS_SHIFT = ERRORS_BITS;
    static constexpr unsigned STRAND_SHIFT = CONTIG_SIZE_BITS + ERRORS_BITS;
    static constexpr unsigned CONTIG_ID_SHIFT = 1 + CONTIG_SIZE_BITS + ERRORS_BITS;
};

template <typename TSpec>
inline uint64_t getSortKey(Match<TSpec> const & me, ContigBegin)
{
    typedef SortKeyShifts<TSpec> Shifts;
    return ((uint64_t)getMember(me, ContigId())    << Shifts::CONTIG_ID_SHIFT) |
           ((uint64_t)onReverseStrand(me)          << Shifts::STRAND_SHIFT)    |
           ((uint64_t)getMember(me, ContigBegin()) << Shifts::CONTIG_POS_SHIFT)|
           ((uint64_t)getMember(me, Errors()));
}
```

**Rationale:** Compiler should optimize this anyway, but making it explicit ensures no runtime computation and improves code readability.

**Risk:** Very low - compile-time only change.

---

### 1.3 Reserve Match Vector Capacity

**File:** `mapper.h:764`

**Current:**
```cpp
reserve(me.matchesByCoord, countHits(me) / 3);
```

**Enhancement:** Add capacity tracking to avoid repeated reallocation:
```cpp
// Before processing loop, estimate capacity more accurately
size_t estimatedMatches = std::max(
    countHits(me) / 2,  // More conservative estimate
    length(me.matchesByCoord)  // Don't shrink if already larger
);
reserve(me.matchesByCoord, estimatedMatches);
```

**Risk:** Low - may use slightly more memory.

---

## Phase 2: Threading Improvements (Weeks 2-3)

### 2.1 Thread-Local Match Buffers

**Files:** `mapper_extender.h`, `mapper.h`

**Problem:** `appendValue()` at line 304 uses atomic operations for every match in parallel mode.

**Solution:** Implement thread-local accumulation with periodic bulk flush.

**Changes to `mapper_extender.h`:**

```cpp
template <typename TSpec, typename Traits>
struct HitsExtender
{
    // ... existing members ...

    // NEW: Thread-local match buffer
    typedef typename Traits::TMatches TLocalMatches;
    TLocalMatches localMatches;
    static constexpr size_t LOCAL_BUFFER_SIZE = 1024;

    // In constructor, reserve local buffer
    HitsExtender(...) : ... {
        reserve(localMatches, LOCAL_BUFFER_SIZE);
        _extendReadsImpl(*this);
        _flushLocalMatches(*this);  // Final flush
    }
};

template <typename TSpec, typename Traits, typename TMatchPos, typename TMatchErrors>
inline void _addMatchImpl(HitsExtender<TSpec, Traits> & me, ...)
{
    setContigPosition(me.prototype, matchBegin, matchEnd);
    me.prototype.errors = matchErrors;

    // Append to local buffer
    appendValue(me.localMatches, me.prototype);

    // Flush when buffer is full
    if (length(me.localMatches) >= HitsExtender<TSpec, Traits>::LOCAL_BUFFER_SIZE)
        _flushLocalMatches(me);

    // ... rest unchanged ...
}

template <typename TSpec, typename Traits>
inline void _flushLocalMatches(HitsExtender<TSpec, Traits> & me)
{
    if (empty(me.localMatches)) return;

    SEQAN_OMP_PRAGMA(critical(HitsExtender_flush))
    {
        append(me.matches, me.localMatches);
    }
    clear(me.localMatches);
}
```

**Risk:** Medium - changes data flow, requires careful testing.

**Validation:**
- Verify match counts are identical
- Verify all matches present (order may differ before sorting)

---

### 2.2 Parallelize Bucket Processing

**File:** `mapper.h:689-690`

**Current:**
```cpp
for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
    TSeedsRanker ranker(hitsCounts, me.ranks[bucketId], ...);
```

**Change:**
```cpp
// Each bucket is independent, can process in parallel
SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
{
    typename TTraits::THitsCounts localHitsCounts;
    TSeedsRanker ranker(localHitsCounts, me.ranks[bucketId], ...);
}
```

**Note:** `hitsCounts` appears to be used only locally within each ranker, so this should be safe. Verify by reading `mapper_ranker.h`.

**Risk:** Medium - need to verify no shared state between buckets.

---

### 2.3 Lock-Free Output Queue

**File:** `mapper_writer.h`

**Problem:** Critical section at line 671 serializes all output.

**Solution:** Replace critical section with lock-free queue and dedicated writer thread.

**New architecture:**

```cpp
// New class for async writing
template <typename TOutputFile>
struct AsyncWriter
{
    TOutputFile & outputFile;
    std::atomic<bool> done{false};

    // Lock-free queue (use boost::lockfree::queue or similar)
    // Or use double-buffering approach
    std::vector<CharString> buffers[2];
    std::atomic<int> activeBuffer{0};
    std::mutex flushMutex;
    std::condition_variable flushCV;
    std::thread writerThread;

    AsyncWriter(TOutputFile & file) : outputFile(file) {
        writerThread = std::thread([this]() { writerLoop(); });
    }

    void submit(CharString && buffer) {
        int idx = activeBuffer.load();
        {
            std::lock_guard<std::mutex> lock(flushMutex);
            buffers[idx].push_back(std::move(buffer));
        }
        flushCV.notify_one();
    }

    void writerLoop() {
        while (!done.load() || !buffersEmpty()) {
            std::unique_lock<std::mutex> lock(flushMutex);
            flushCV.wait(lock, [this]{ return done.load() || !buffersEmpty(); });

            int idx = activeBuffer.load();
            activeBuffer.store(1 - idx);  // Swap buffers
            lock.unlock();

            // Write old buffer
            for (auto & buf : buffers[1-idx]) {
                write(outputFile.stream, buf);
            }
            buffers[1-idx].clear();
        }
    }

    ~AsyncWriter() {
        done.store(true);
        flushCV.notify_one();
        writerThread.join();
    }
};
```

**Risk:** High - significant architectural change.

**Alternative (lower risk):** Simply increase buffer size (Phase 1.1) and accept some serialization.

---

## Phase 3: Algorithm Optimizations (Weeks 3-4)

### 3.1 Early Termination in Extension Loop

**File:** `mapper_extender.h:263-285`

**Current:** Loop continues through all SA positions even if read is already well-mapped.

**Enhancement:** Add early exit when sufficient matches found:

```cpp
for (TSAPos saPos = getValueI1(hitRange); saPos < getValueI2(hitRange); ++saPos)
{
    // Early termination: if we have enough matches at optimal error level, stop
    TReadId readId = getReadId(me.readSeqs, readSeqId);
    if (isMapped(me.ctx, readId)) {
        // Check if we have enough matches in best stratum
        // This requires tracking match count per read
        break;
    }

    // ... existing extension code ...
}
```

**Risk:** Medium - must not break correctness for strata-count reporting.

**Validation:** Compare X0/X1 tags in output.

---

### 3.2 Batch SA Value Prefetching

**File:** `mapper_extender.h:263-285`

**Optimization:** Prefetch SA values to hide memory latency:

```cpp
// Before loop, prefetch first few SA values
constexpr size_t PREFETCH_DISTANCE = 8;
for (TSAPos i = getValueI1(hitRange);
     i < std::min(getValueI1(hitRange) + PREFETCH_DISTANCE, getValueI2(hitRange));
     ++i) {
    __builtin_prefetch(&me.sa[i], 0, 1);
}

for (TSAPos saPos = getValueI1(hitRange); saPos < getValueI2(hitRange); ++saPos)
{
    // Prefetch ahead
    if (saPos + PREFETCH_DISTANCE < getValueI2(hitRange))
        __builtin_prefetch(&me.sa[saPos + PREFETCH_DISTANCE], 0, 1);

    TSAValue saValue = me.sa[saPos];
    // ... rest unchanged ...
}
```

**Risk:** Low - CPU hint only, no functional change.

**Note:** Effectiveness depends on hardware. May need `#ifdef` for portability.

---

### 3.3 Optimize LCP Computation in Extension

**File:** `find_extender.h:129-136, 195-202`

**Current:** LCP computed with temporary reversed strings:
```cpp
TNeedleInfixRev needleInfixRev(needleInfix);
THaystackInfixRev haystackInfixRev(haystackInfix);
lcp = lcpLength(haystackInfixRev, needleInfixRev);
```

**Optimization:** Implement reverse LCP without constructing reversed strings:

```cpp
template <typename TSeq1, typename TSeq2>
inline typename Size<TSeq1>::Type lcpLengthReverse(TSeq1 const & seq1, TSeq2 const & seq2)
{
    typedef typename Size<TSeq1>::Type TSize;
    TSize len = std::min(length(seq1), length(seq2));
    TSize i = 0;
    while (i < len && seq1[length(seq1) - 1 - i] == seq2[length(seq2) - 1 - i])
        ++i;
    return i;
}
```

**Risk:** Low - equivalent computation.

---

## Phase 4: Memory Layout Optimizations (Week 4-5)

### 4.1 Match Structure Unpacking for Hot Paths

**File:** `bits_matches.h`

**Problem:** Packed bitfields require byte-level access, poor for CPU pipelines.

**Solution:** Create unpacked view for sorting/filtering operations:

```cpp
// Unpacked match for fast comparison
template <typename TSpec>
struct MatchUnpacked
{
    uint32_t readId;
    uint32_t contigId;
    uint32_t contigBegin;
    uint32_t contigEnd;
    uint16_t errors;
    uint8_t flags;  // strand, validity

    MatchUnpacked() = default;
    explicit MatchUnpacked(Match<TSpec> const & m) {
        readId = getMember(m, ReadId());
        contigId = getMember(m, ContigId());
        contigBegin = getMember(m, ContigBegin());
        contigEnd = getMember(m, ContigEnd());
        errors = getMember(m, Errors());
        flags = onReverseStrand(m) ? 1 : 0;
        if (isValid(m)) flags |= 2;
    }
};

// Use for sorting, then pack back
template <typename TMatches>
void sortMatchesFast(TMatches & matches) {
    std::vector<MatchUnpacked<...>> unpacked;
    unpacked.reserve(length(matches));
    for (auto const & m : matches)
        unpacked.emplace_back(m);

    std::sort(unpacked.begin(), unpacked.end(), ...);

    // Pack back
    for (size_t i = 0; i < unpacked.size(); ++i)
        packMatch(matches[i], unpacked[i]);
}
```

**Risk:** Medium - requires careful bit preservation.

**Trade-off:** Uses 2x memory during sort. Only worthwhile for large match sets.

---

### 4.2 Hit Structure Reordering

**File:** `bits_hits.h:59-65`

**Current:**
```cpp
struct Hit<TSize, HammingDistance> {
    typename Position<Hit>::Type range;     // Pair<uint32_t> = 8 bytes
    typename Id<Hit>::Type seedId;          // unsigned = 4 bytes
    unsigned char errors;                    // 1 byte + 3 padding
};  // Total: 16 bytes due to alignment
```

**Optimized:**
```cpp
struct Hit<TSize, HammingDistance> {
    typename Position<Hit>::Type range;     // 8 bytes
    typename Id<Hit>::Type seedId;          // 4 bytes
    unsigned char errors;                    // 1 byte
    unsigned char padding[3];               // explicit padding
};  // Same size, but explicit
```

**Alternative (if seedId can be smaller):**
```cpp
struct Hit<TSize, HammingDistance> {
    typename Position<Hit>::Type range;     // 8 bytes
    uint16_t seedId;                        // 2 bytes (if < 65k seeds)
    uint8_t errors;                         // 1 byte
    uint8_t reserved;                       // 1 byte
};  // Total: 12 bytes
```

**Risk:** Low for explicit padding, Medium for size reduction.

---

## Phase 5: I/O Optimizations (Week 5-6)

### 5.1 Parallel File Prefetching

**File:** `file_prefetched.h`

**Current:** Single reader thread.

**Enhancement:** Multi-block prefetching with larger buffers:

```cpp
template <typename TFile, typename TRecords, typename TSize>
struct PrefetchedFile<TFile, TRecords, TSize, Parallel>
{
    // ... existing members ...

    static constexpr size_t NUM_BUFFERS = 3;  // Triple buffering
    TRecords buffers[NUM_BUFFERS];
    std::atomic<int> readIdx{0};
    std::atomic<int> writeIdx{0};
    std::atomic<int> processIdx{0};

    // Reader thread fills buffers ahead
    void readerLoop() {
        while (!eof) {
            int idx = writeIdx.load();
            if ((idx + 1) % NUM_BUFFERS == processIdx.load())
                continue;  // Buffer full, wait

            readRecords(buffers[idx], file, maxRecords);
            writeIdx.store((idx + 1) % NUM_BUFFERS);
        }
    }
};
```

**Risk:** Medium - threading complexity.

---

### 5.2 Memory-Mapped Index Loading

**File:** `mapper.h` (index loading)

**If not already using mmap:** Ensure index files are memory-mapped rather than read into RAM.

**Check:** Verify `MMap<>` template parameter is being used effectively.

---

## Phase 6: Build System Optimizations

### 6.1 Profile-Guided Optimization (PGO)

**File:** `CMakeLists.txt`

**Add PGO support:**
```cmake
option(YARA_PGO_GENERATE "Build for PGO profile generation" OFF)
option(YARA_PGO_USE "Build using PGO profile" OFF)

if (YARA_PGO_GENERATE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate")
endif()

if (YARA_PGO_USE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-use")
endif()
```

**Usage:**
```bash
# Generate profile
cmake -DYARA_PGO_GENERATE=ON ..
make
./yara_mapper ... # Run typical workload
# Use profile
cmake -DYARA_PGO_USE=ON ..
make
```

---

### 6.2 Link-Time Optimization (LTO)

**File:** `CMakeLists.txt`

```cmake
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
    set_property(TARGET yara_mapper PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET yara_indexer PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
```

---

## Validation Plan

### Functional Testing
1. Run `python tests/run_tests.py` after each change
2. Compare SAM output byte-for-byte (excluding @PG line)
3. Verify X0, X1, NM tags are identical

### Performance Testing
1. Create benchmark script with timing
2. Test with various thread counts (1, 2, 4, 8, 16)
3. Test with HLA reference + typical read sets
4. Measure: wall time, CPU time, peak memory, I/O wait

### Regression Testing
1. Test on multiple platforms (Linux, macOS)
2. Test with different compilers (GCC, Clang)
3. Test single-end and paired-end modes

---

## Implementation Order (Recommended)

1. **Week 1:** Phases 1.1, 1.2, 1.3, 6.2 (quick wins + LTO)
2. **Week 2:** Phase 2.1 (thread-local buffers)
3. **Week 3:** Phases 2.2, 3.1, 3.2 (parallelization + early termination)
4. **Week 4:** Phases 3.3, 4.2 (algorithm + memory)
5. **Week 5:** Phase 5.1 or 2.3 (I/O - choose based on profiling)
6. **Week 6:** Phase 6.1 (PGO), final testing, documentation

---

## Risk Summary

| Phase | Risk | Mitigation |
|-------|------|------------|
| 1.x | Very Low | Simple constant changes |
| 2.1 | Medium | Extensive testing of match completeness |
| 2.2 | Medium | Verify bucket independence |
| 2.3 | High | Consider as optional; buffer increase may suffice |
| 3.x | Low-Medium | Functional tests catch regressions |
| 4.x | Medium | Bit-level validation |
| 5.x | Medium | I/O correctness tests |
| 6.x | Low | Build-only changes |
