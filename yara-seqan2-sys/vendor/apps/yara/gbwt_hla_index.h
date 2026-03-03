// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2026, Enrico Siragusa, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// GBWT-based HLA index for finding and extension optimization
// ==========================================================================

#ifndef APP_YARA_GBWT_HLA_INDEX_H_
#define APP_YARA_GBWT_HLA_INDEX_H_

#ifdef YARA_WITH_GBWT

#include <gbwt/gbwt.h>
#include <gbwt/dynamic_gbwt.h>
#include <sdsl/io.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/index.h>

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <mutex>

using namespace seqan2;

namespace yara_gbwt {

// ============================================================================
// FMIndex Configuration for HLA sequences
// ============================================================================
// Simple FMIndex configuration for the HLA allele sequence index

struct HLAFMConfig
{
    typedef HLAFMConfig                                     TMe;

    // Text configuration
    typedef Owner<ConcatDirect<TMe> >                       TSSetSpec_;
    typedef StringSet<String<Dna5, Alloc<> >, TSSetSpec_>   Text;

    // Length sum type (use 64-bit for potentially large allele collections)
    typedef uint64_t                                        LengthSum;

    // RankDictionary configuration
    typedef Levels<void, TMe>                               Bwt;
    typedef Levels<void, TMe>                               Sentinels;

    // RankDictionary requirements
    typedef Alloc<>                                         Fibre;
    typedef uint64_t                                        Size;

    // Sparse SA sampling configuration
    static const unsigned SAMPLING =                        10;
    static const unsigned WORDS_PER_BLOCK =                 1;
    static const unsigned LEVELS =                          1;
};

} // namespace yara_gbwt

namespace seqan2 {
// SAValue specialization for HLA sequence index
template <typename TValue, typename TSpec>
struct SAValue<StringSet<String<TValue, TSpec>, Owner<ConcatDirect<yara_gbwt::HLAFMConfig> > > >
{
    typedef Pair<uint32_t, uint32_t, Pack>   Type;
};

// _getNodeByChar overloads for HLAFMConfig
// These are required for FM-index search operations

template <typename TText, typename TOccSpec, typename TSpec, typename TSize>
inline bool
_getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> >, VSTree<TopDown<TSpec> > > const & it,
               typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> > >::Type const & vDesc,
               Pair<typename Size<Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> > >::Type> & _range,
               TSize & /*smaller*/,
               Dna5 c)
{
    typedef Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> >  TIndex;
    typedef typename Fibre<TIndex, FibreLF>::Type                    TLF;
    typedef typename Value<TIndex>::Type                             TAlphabet;

    TIndex const & index = container(it);
    TLF const & lf = indexLF(index);

    if (ordValue(c) >= ValueSize<TAlphabet>::VALUE) return false;

    _range = range(index, vDesc);
    _range.i1 = lf(_range.i1, c);
    _range.i2 = lf(_range.i2, c);

    return _range.i1 < _range.i2;
}

template <typename TText, typename TOccSpec, typename TSpec, typename TSize>
inline bool
_getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> >, VSTree<TopDown<TSpec> > > const & it,
               typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> > >::Type const & vDesc,
               Pair<typename Size<Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> > >::Type> & _range,
               TSize & /*smaller*/,
               Dna5Q c)
{
    typedef Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> >  TIndex;
    typedef typename Fibre<TIndex, FibreLF>::Type                    TLF;
    typedef typename Value<TIndex>::Type                             TAlphabet;

    TIndex const & index = container(it);
    TLF const & lf = indexLF(index);

    if (ordValue(c) >= ValueSize<TAlphabet>::VALUE) return false;

    _range = range(index, vDesc);
    _range.i1 = lf(_range.i1, c);
    _range.i2 = lf(_range.i2, c);

    return _range.i1 < _range.i2;
}

// Generic _getNodeByChar for any character type
template <typename TText, typename TOccSpec, typename TSpec, typename TSize, typename TChar>
inline bool
_getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> >, VSTree<TopDown<TSpec> > > const & it,
               typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> > >::Type const & vDesc,
               Pair<typename Size<Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> > >::Type> & _range,
               TSize & /*smaller*/,
               TChar c)
{
    typedef Index<TText, FMIndex<TOccSpec, yara_gbwt::HLAFMConfig> >  TIndex;
    typedef typename Fibre<TIndex, FibreLF>::Type                    TLF;

    TIndex const & index = container(it);
    TLF const & lf = indexLF(index);

    _range = range(index, vDesc);
    _range.i1 = lf(_range.i1, c);
    _range.i2 = lf(_range.i2, c);

    return _range.i1 < _range.i2;
}

} // namespace seqan2

namespace yara_gbwt {

// ============================================================================
// Constants
// ============================================================================

static constexpr size_t MIN_CONSERVED_LENGTH = 5;
static constexpr size_t GBWT_BATCH_SIZE = 1000;
static constexpr double REGION_DETECTION_IDENTITY = 0.90;  // 90% identity threshold

// ============================================================================
// HLARegionMap - Maps reference positions to HLA allele graph
// ============================================================================

struct HLARegionMap
{
    struct Region {
        uint32_t contigId;
        uint32_t startPos;
        uint32_t endPos;
        uint32_t alleleSetId;
        CharString geneName;
    };

    std::vector<Region> regions;
    std::map<uint32_t, std::vector<size_t>> contigRegions;

    bool isHLARegion(uint32_t contigId, uint32_t pos) const
    {
        auto it = contigRegions.find(contigId);
        if (it == contigRegions.end()) return false;
        for (size_t idx : it->second)
        {
            if (pos >= regions[idx].startPos && pos < regions[idx].endPos)
                return true;
        }
        return false;
    }

    Region const * getRegion(uint32_t contigId, uint32_t pos) const
    {
        auto it = contigRegions.find(contigId);
        if (it == contigRegions.end()) return nullptr;
        for (size_t idx : it->second)
        {
            if (pos >= regions[idx].startPos && pos < regions[idx].endPos)
                return &regions[idx];
        }
        return nullptr;
    }

    void addRegion(uint32_t contigId, uint32_t start, uint32_t end,
                   uint32_t alleleSetId, CharString const & gene)
    {
        size_t idx = regions.size();
        regions.push_back({contigId, start, end, alleleSetId, gene});
        contigRegions[contigId].push_back(idx);
    }

    void sortRegions()
    {
        for (auto & pair : contigRegions)
        {
            std::sort(pair.second.begin(), pair.second.end(),
                [this](size_t a, size_t b) {
                    return regions[a].startPos < regions[b].startPos;
                });
        }
    }

    void clear()
    {
        regions.clear();
        contigRegions.clear();
    }

    bool save(const char * fileName) const;
    bool load(const char * fileName);
};

// ============================================================================
// AlleleCoordMap - Maps positions in allele sequences to reference coordinates
// ============================================================================

struct AlleleCoordMap
{
    // For each allele, store (refContigId, refStart, length, strand)
    struct AlleleInfo {
        uint32_t alleleId;
        uint32_t alleleSetId;
        uint32_t refContigId;
        uint32_t refStart;
        uint32_t refEnd;
        bool reverseStrand;
    };

    std::vector<AlleleInfo> alleles;

    // Index: alleleSetId -> list of allele indices for efficient lookup
    std::map<uint32_t, std::vector<size_t>> allelesBySet;

    // Map from (alleleSetId, posInAllele) -> (refContigId, refPos)
    // For indels, we need per-base coordinate mapping
    // Simplified: store just the start and assume linear mapping
    // For full indel support, would need CIGAR-based mapping

    void addAllele(uint32_t alleleId, uint32_t alleleSetId,
                   uint32_t refContigId, uint32_t refStart, uint32_t refEnd,
                   bool reverseStrand)
    {
        size_t idx = alleles.size();
        alleles.push_back({alleleId, alleleSetId, refContigId, refStart, refEnd, reverseStrand});
        allelesBySet[alleleSetId].push_back(idx);
    }

    // Get all alleles for a given alleleSetId (efficient O(1) lookup)
    std::vector<size_t> const * getAllelesForSet(uint32_t alleleSetId) const
    {
        auto it = allelesBySet.find(alleleSetId);
        return (it != allelesBySet.end()) ? &it->second : nullptr;
    }

    // Rebuild the index after loading
    void rebuildIndex()
    {
        allelesBySet.clear();
        for (size_t i = 0; i < alleles.size(); ++i)
        {
            allelesBySet[alleles[i].alleleSetId].push_back(i);
        }
    }

    // Translate allele position to reference position
    // Returns false if translation not possible (e.g., indel region)
    bool toRefCoord(uint32_t alleleId, uint32_t posInAllele,
                    uint32_t & refContigId, uint32_t & refPos) const
    {
        if (alleleId >= alleles.size()) return false;
        auto const & info = alleles[alleleId];
        refContigId = info.refContigId;

        if (info.reverseStrand)
        {
            // Reverse strand: position counts from end
            if (posInAllele > info.refEnd - info.refStart) return false;
            refPos = info.refEnd - posInAllele;
        }
        else
        {
            refPos = info.refStart + posInAllele;
            if (refPos > info.refEnd) return false;
        }
        return true;
    }

    bool save(const char * fileName) const;
    bool load(const char * fileName);
};

// ============================================================================
// HLASequenceIndex - DNA sequence index for HLA alleles
// ============================================================================
// Uses a concatenated sequence with delimiter to enable FM-index search
// across all alleles simultaneously.

struct HLASequenceIndex
{
    typedef StringSet<String<Dna5>, Owner<ConcatDirect<HLAFMConfig>>> TAlleleSeqs;
    typedef Index<TAlleleSeqs, FMIndex<void, HLAFMConfig>> TIndex;

    TAlleleSeqs alleleSeqs;
    TIndex index;
    std::vector<uint32_t> alleleSetIds;  // alleleSetId for each allele

    bool build()
    {
        if (empty(alleleSeqs)) return false;
        index = TIndex(alleleSeqs);
        indexRequire(index, FibreSA());
        return true;
    }

    void addAllele(String<Dna5> const & seq, uint32_t alleleSetId)
    {
        appendValue(alleleSeqs, seq);
        alleleSetIds.push_back(alleleSetId);
    }

    size_t alleleCount() const { return length(alleleSeqs); }

    bool save(const char * fileName) const;
    bool load(const char * fileName);
};

// ============================================================================
// HLAGraphBuilder - Constructs GBWT and sequence index from HLA alleles
// ============================================================================

class HLAGraphBuilder
{
public:
    bool loadFromFasta(CharString const & filename, bool verbose = false);

    // Auto-detect HLA regions by aligning first allele of each gene to reference
    template <typename TContigSeqs>
    bool detectReferenceRegions(TContigSeqs const & contigs,
                                StringSet<CharString> const & contigNames,
                                bool verbose = false);

    bool buildIndex(bool verbose = false);
    bool save(CharString const & prefix) const;

    // Load pre-built index
    bool load(CharString const & prefix);

    // Accessors
    gbwt::GBWT const & getGBWT() const { return gbwt_; }
    HLARegionMap const & getRegionMap() const { return regionMap_; }
    AlleleCoordMap const & getCoordMap() const { return coordMap_; }
    HLASequenceIndex const & getSeqIndex() const { return seqIndex_; }

    void setReferenceRegion(uint32_t contigId, uint32_t start, uint32_t end,
                            CharString const & geneName);

private:
    struct AlleleSet {
        CharString geneName;
        StringSet<String<Dna5>> sequences;
        StringSet<CharString> names;
        uint32_t refContigId = 0;
        uint32_t refStart = 0;
        uint32_t refEnd = 0;
    };

    std::vector<AlleleSet> alleleSets_;
    gbwt::GBWT gbwt_;
    HLARegionMap regionMap_;
    AlleleCoordMap coordMap_;
    HLASequenceIndex seqIndex_;

    void buildGBWTPaths_(bool verbose);
};

// ============================================================================
// Implementation: HLARegionMap serialization
// ============================================================================

inline bool HLARegionMap::save(const char * fileName) const
{
    std::ofstream out(fileName, std::ios::binary);
    if (!out) return false;

    size_t count = regions.size();
    out.write(reinterpret_cast<const char*>(&count), sizeof(count));

    for (auto const & r : regions)
    {
        out.write(reinterpret_cast<const char*>(&r.contigId), sizeof(r.contigId));
        out.write(reinterpret_cast<const char*>(&r.startPos), sizeof(r.startPos));
        out.write(reinterpret_cast<const char*>(&r.endPos), sizeof(r.endPos));
        out.write(reinterpret_cast<const char*>(&r.alleleSetId), sizeof(r.alleleSetId));

        size_t nameLen = length(r.geneName);
        out.write(reinterpret_cast<const char*>(&nameLen), sizeof(nameLen));
        out.write(toCString(r.geneName), nameLen);
    }
    return out.good();
}

inline bool HLARegionMap::load(const char * fileName)
{
    std::ifstream in(fileName, std::ios::binary);
    if (!in) return false;

    clear();

    size_t count;
    in.read(reinterpret_cast<char*>(&count), sizeof(count));

    for (size_t i = 0; i < count; ++i)
    {
        Region r;
        in.read(reinterpret_cast<char*>(&r.contigId), sizeof(r.contigId));
        in.read(reinterpret_cast<char*>(&r.startPos), sizeof(r.startPos));
        in.read(reinterpret_cast<char*>(&r.endPos), sizeof(r.endPos));
        in.read(reinterpret_cast<char*>(&r.alleleSetId), sizeof(r.alleleSetId));

        size_t nameLen;
        in.read(reinterpret_cast<char*>(&nameLen), sizeof(nameLen));
        resize(r.geneName, nameLen);
        in.read(&r.geneName[0], nameLen);

        size_t idx = regions.size();
        regions.push_back(r);
        contigRegions[r.contigId].push_back(idx);
    }

    sortRegions();
    return in.good();
}

// ============================================================================
// Implementation: AlleleCoordMap serialization
// ============================================================================

inline bool AlleleCoordMap::save(const char * fileName) const
{
    std::ofstream out(fileName, std::ios::binary);
    if (!out) return false;

    size_t count = alleles.size();
    out.write(reinterpret_cast<const char*>(&count), sizeof(count));

    for (auto const & a : alleles)
    {
        out.write(reinterpret_cast<const char*>(&a.alleleId), sizeof(a.alleleId));
        out.write(reinterpret_cast<const char*>(&a.alleleSetId), sizeof(a.alleleSetId));
        out.write(reinterpret_cast<const char*>(&a.refContigId), sizeof(a.refContigId));
        out.write(reinterpret_cast<const char*>(&a.refStart), sizeof(a.refStart));
        out.write(reinterpret_cast<const char*>(&a.refEnd), sizeof(a.refEnd));
        uint8_t rev = a.reverseStrand ? 1 : 0;
        out.write(reinterpret_cast<const char*>(&rev), sizeof(rev));
    }
    return out.good();
}

inline bool AlleleCoordMap::load(const char * fileName)
{
    std::ifstream in(fileName, std::ios::binary);
    if (!in) return false;

    alleles.clear();
    allelesBySet.clear();

    size_t count;
    in.read(reinterpret_cast<char*>(&count), sizeof(count));

    for (size_t i = 0; i < count; ++i)
    {
        AlleleInfo a;
        in.read(reinterpret_cast<char*>(&a.alleleId), sizeof(a.alleleId));
        in.read(reinterpret_cast<char*>(&a.alleleSetId), sizeof(a.alleleSetId));
        in.read(reinterpret_cast<char*>(&a.refContigId), sizeof(a.refContigId));
        in.read(reinterpret_cast<char*>(&a.refStart), sizeof(a.refStart));
        in.read(reinterpret_cast<char*>(&a.refEnd), sizeof(a.refEnd));
        uint8_t rev;
        in.read(reinterpret_cast<char*>(&rev), sizeof(rev));
        a.reverseStrand = (rev != 0);
        alleles.push_back(a);
    }

    // Rebuild the index for efficient lookup
    rebuildIndex();

    return in.good();
}

// ============================================================================
// Implementation: HLASequenceIndex serialization
// ============================================================================

inline bool HLASequenceIndex::save(const char * fileName) const
{
    // Save alleleSetIds to main file
    std::ofstream out(fileName, std::ios::binary);
    if (!out) return false;

    size_t count = alleleSetIds.size();
    out.write(reinterpret_cast<const char*>(&count), sizeof(count));

    for (size_t i = 0; i < count; ++i)
    {
        out.write(reinterpret_cast<const char*>(&alleleSetIds[i]), sizeof(alleleSetIds[i]));
    }

    out.close();

    // Save sequences (text) - needed for coordinate lookups
    String<char> textFile = fileName;
    append(textFile, ".txt");
    if (!seqan2::save(alleleSeqs, toCString(textFile)))
        return false;

    // Save FM-index components
    String<char> saFile = fileName;
    append(saFile, ".sa");
    if (!seqan2::save(getFibre(index, FibreSA()), toCString(saFile)))
        return false;

    String<char> lfFile = fileName;
    append(lfFile, ".lf");
    if (!seqan2::save(getFibre(index, FibreLF()), toCString(lfFile)))
        return false;

    return true;
}

inline bool HLASequenceIndex::load(const char * fileName)
{
    // Load alleleSetIds from main file
    std::ifstream in(fileName, std::ios::binary);
    if (!in) return false;

    size_t count;
    in.read(reinterpret_cast<char*>(&count), sizeof(count));

    alleleSetIds.resize(count);
    for (size_t i = 0; i < count; ++i)
    {
        in.read(reinterpret_cast<char*>(&alleleSetIds[i]), sizeof(alleleSetIds[i]));
    }
    in.close();

    // Load sequences (text)
    String<char> textFile = fileName;
    append(textFile, ".txt");
    if (!seqan2::open(alleleSeqs, toCString(textFile)))
        return false;

    // Rebuild FM-index from loaded sequences
    // Note: We can't just load the SA and LF fibres separately because the
    // FM-index needs to be constructed with a reference to the text.
    // Rebuilding is the safest approach.
    if (empty(alleleSeqs))
        return false;

    index = TIndex(alleleSeqs);
    indexRequire(index, FibreSA());

    return true;
}

// ============================================================================
// Implementation: HLAGraphBuilder
// ============================================================================

inline bool HLAGraphBuilder::loadFromFasta(CharString const & filename, bool verbose)
{
    SeqFileIn seqFile;
    if (!open(seqFile, toCString(filename)))
    {
        std::cerr << "Error: Cannot open HLA allele file: " << filename << std::endl;
        return false;
    }

    std::map<CharString, size_t> geneToIndex;

    CharString id;
    String<Dna5> seq;

    while (!atEnd(seqFile))
    {
        readRecord(id, seq, seqFile);

        // Parse gene name from ID (e.g., "HLA-A*01:01" -> "HLA-A")
        CharString geneName;
        for (size_t i = 0; i < length(id) && id[i] != '*'; ++i)
            appendValue(geneName, id[i]);

        auto it = geneToIndex.find(geneName);
        if (it == geneToIndex.end())
        {
            geneToIndex[geneName] = alleleSets_.size();
            alleleSets_.push_back({geneName, {}, {}, 0, 0, 0});
        }
        size_t idx = geneToIndex[geneName];

        appendValue(alleleSets_[idx].sequences, seq);
        appendValue(alleleSets_[idx].names, id);
    }

    if (verbose)
    {
        std::cerr << "Loaded " << alleleSets_.size() << " HLA genes:" << std::endl;
        for (auto const & as : alleleSets_)
            std::cerr << "  " << as.geneName << ": " << length(as.sequences) << " alleles" << std::endl;
    }

    return !alleleSets_.empty();
}

template <typename TContigSeqs>
inline bool HLAGraphBuilder::detectReferenceRegions(
    TContigSeqs const & contigs,
    StringSet<CharString> const & contigNames,
    bool verbose)
{
    if (alleleSets_.empty()) return false;

    // Build k-mer index for fast candidate position finding
    static constexpr size_t KMER_SIZE = 15;
    std::unordered_map<uint64_t, std::vector<std::pair<uint32_t, uint32_t>>> kmerIndex;

    if (verbose)
        std::cerr << "  Building k-mer index for reference contigs..." << std::endl;

    // Hash function for k-mers
    auto hashKmer = [](const auto & seq, size_t pos, size_t k) -> uint64_t {
        uint64_t hash = 0;
        for (size_t i = 0; i < k; ++i)
        {
            hash = hash * 5 + ordValue(seq[pos + i]);
        }
        return hash;
    };

    // Index all contigs
    for (size_t contigId = 0; contigId < length(contigs); ++contigId)
    {
        auto const & contig = contigs[contigId];
        if (length(contig) < KMER_SIZE) continue;

        for (size_t pos = 0; pos + KMER_SIZE <= length(contig); pos += 50)  // Sample every 50bp
        {
            uint64_t hash = hashKmer(contig, pos, KMER_SIZE);
            kmerIndex[hash].emplace_back(contigId, pos);
        }
    }

    if (verbose)
        std::cerr << "  Searching for " << alleleSets_.size() << " HLA genes in parallel..." << std::endl;

    std::mutex outputMutex;

    #pragma omp parallel for schedule(dynamic)
    for (size_t asIdx = 0; asIdx < alleleSets_.size(); ++asIdx)
    {
        auto & as = alleleSets_[asIdx];
        if (empty(as.sequences)) continue;

        // Use the first (typically reference) allele for detection
        auto const & refAllele = as.sequences[0];
        size_t refAlleleLen = length(refAllele);
        if (refAlleleLen < KMER_SIZE) continue;

        uint32_t bestContigId = 0;
        uint32_t bestStart = 0;
        uint32_t bestEnd = 0;
        size_t bestScore = 0;

        // Get candidate positions from k-mer index
        std::set<std::pair<uint32_t, uint32_t>> candidates;

        // Sample k-mers from the allele
        for (size_t pos = 0; pos + KMER_SIZE <= refAlleleLen; pos += refAlleleLen / 5)
        {
            uint64_t hash = hashKmer(refAllele, pos, KMER_SIZE);
            auto it = kmerIndex.find(hash);
            if (it != kmerIndex.end())
            {
                for (auto const & loc : it->second)
                {
                    // Adjust position to potential allele start
                    int32_t adjustedStart = static_cast<int32_t>(loc.second) - static_cast<int32_t>(pos);
                    if (adjustedStart >= 0)
                    {
                        candidates.insert({loc.first, static_cast<uint32_t>(adjustedStart)});
                    }
                }
            }
        }

        // Verify candidates
        for (auto const & cand : candidates)
        {
            uint32_t contigId = cand.first;
            uint32_t startPos = cand.second;

            auto const & contig = contigs[contigId];
            if (startPos + refAlleleLen > length(contig)) continue;

            // Count matches
            size_t matches = 0;
            for (size_t i = 0; i < refAlleleLen; ++i)
            {
                if (contig[startPos + i] == refAllele[i])
                    ++matches;
            }

            if (matches > bestScore)
            {
                bestScore = matches;
                bestContigId = contigId;
                bestStart = startPos;
                bestEnd = startPos + refAlleleLen;

                // Early termination for near-perfect match
                if (static_cast<double>(matches) / refAlleleLen >= 0.99)
                    break;
            }
        }

        // If k-mer search didn't find good match, fall back to gene-name based search
        if (bestScore == 0 || static_cast<double>(bestScore) / refAlleleLen < REGION_DETECTION_IDENTITY)
        {
            // Try to find contig with matching gene name
            for (size_t contigId = 0; contigId < length(contigs); ++contigId)
            {
                CharString const & name = contigNames[contigId];
                // Check if contig name contains gene name (e.g., "HLA-A" in "HLA00001")
                if (length(name) >= 3 && length(as.geneName) >= 5)
                {
                    // Extract gene letter (e.g., 'A' from 'HLA-A')
                    char geneLetter = as.geneName[4];
                    // Check if this is an HLA contig for the same gene
                    bool matchesGene = false;
                    for (size_t i = 0; i + 1 < length(name); ++i)
                    {
                        if (name[i] == 'H' && name[i+1] == 'L' && name[i+2] == 'A')
                        {
                            matchesGene = true;
                            break;
                        }
                    }

                    if (matchesGene)
                    {
                        auto const & contig = contigs[contigId];
                        if (length(contig) < refAlleleLen) continue;

                        // Align at start
                        size_t matches = 0;
                        for (size_t i = 0; i < refAlleleLen && i < length(contig); ++i)
                        {
                            if (contig[i] == refAllele[i])
                                ++matches;
                        }

                        if (matches > bestScore)
                        {
                            bestScore = matches;
                            bestContigId = contigId;
                            bestStart = 0;
                            bestEnd = refAlleleLen;
                        }
                    }
                }
            }
        }

        if (bestScore > 0 && static_cast<double>(bestScore) / refAlleleLen >= REGION_DETECTION_IDENTITY)
        {
            as.refContigId = bestContigId;
            as.refStart = bestStart;
            as.refEnd = bestEnd;

            if (verbose)
            {
                std::lock_guard<std::mutex> lock(outputMutex);
                std::cerr << "  Found " << as.geneName << " at "
                          << contigNames[bestContigId] << ":"
                          << bestStart << "-" << bestEnd
                          << " (" << (100.0 * bestScore / refAlleleLen) << "% identity)"
                          << std::endl;
            }
        }
        else if (verbose)
        {
            std::lock_guard<std::mutex> lock(outputMutex);
            std::cerr << "  Warning: Could not locate " << as.geneName << " in reference" << std::endl;
        }
    }

    return true;
}

inline void HLAGraphBuilder::setReferenceRegion(
    uint32_t contigId, uint32_t start, uint32_t end, CharString const & geneName)
{
    for (auto & as : alleleSets_)
    {
        if (as.geneName == geneName)
        {
            as.refContigId = contigId;
            as.refStart = start;
            as.refEnd = end;
            return;
        }
    }
}

inline void HLAGraphBuilder::buildGBWTPaths_(bool verbose)
{
    // Build GBWT with one path per allele
    // Each path represents a complete allele sequence
    // Nodes are created for each allele (simplest representation)

    size_t totalAlleles = 0;
    for (auto const & as : alleleSets_)
        totalAlleles += length(as.sequences);

    gbwt::GBWTBuilder builder(totalAlleles + 1, GBWT_BATCH_SIZE);

    gbwt::node_type nodeId = 1;
    for (size_t asIdx = 0; asIdx < alleleSets_.size(); ++asIdx)
    {
        auto const & as = alleleSets_[asIdx];
        for (size_t aIdx = 0; aIdx < length(as.sequences); ++aIdx)
        {
            gbwt::vector_type path;
            path.push_back(gbwt::Node::encode(nodeId, false));
            builder.insert(path, true);
            ++nodeId;
        }
    }

    builder.finish();
    gbwt_ = gbwt::GBWT(builder.index);

    if (verbose)
    {
        std::cerr << "Built GBWT with " << gbwt_.sequences() << " paths" << std::endl;
    }
}

inline bool HLAGraphBuilder::buildIndex(bool verbose)
{
    if (alleleSets_.empty()) return false;

    regionMap_.clear();
    coordMap_.alleles.clear();

    uint32_t alleleId = 0;
    for (size_t asIdx = 0; asIdx < alleleSets_.size(); ++asIdx)
    {
        auto & as = alleleSets_[asIdx];

        // Add region if reference coordinates are set
        if (as.refEnd > as.refStart)
        {
            regionMap_.addRegion(as.refContigId, as.refStart, as.refEnd, asIdx, as.geneName);
        }

        // Add each allele to the sequence index and coordinate map
        for (size_t aIdx = 0; aIdx < length(as.sequences); ++aIdx)
        {
            seqIndex_.addAllele(as.sequences[aIdx], asIdx);
            coordMap_.addAllele(alleleId, asIdx, as.refContigId, as.refStart, as.refEnd, false);
            ++alleleId;
        }
    }

    regionMap_.sortRegions();

    // Build sequence FM-index
    if (!seqIndex_.build())
    {
        if (verbose) std::cerr << "Error building sequence index" << std::endl;
        return false;
    }

    // Build GBWT
    buildGBWTPaths_(verbose);

    if (verbose)
    {
        std::cerr << "Built HLA index: " << seqIndex_.alleleCount() << " alleles, "
                  << regionMap_.regions.size() << " regions" << std::endl;
    }

    return true;
}

inline bool HLAGraphBuilder::save(CharString const & prefix) const
{
    // Save GBWT
    CharString gbwtFile = prefix;
    append(gbwtFile, ".gbwt");
    if (!sdsl::store_to_file(gbwt_, toCString(gbwtFile)))
        return false;

    // Save region map
    CharString regionFile = prefix;
    append(regionFile, ".gbwt.regions");
    if (!regionMap_.save(toCString(regionFile)))
        return false;

    // Save coordinate map
    CharString coordFile = prefix;
    append(coordFile, ".gbwt.coords");
    if (!coordMap_.save(toCString(coordFile)))
        return false;

    // Save sequence index
    CharString seqFile = prefix;
    append(seqFile, ".gbwt.seqs");
    if (!seqIndex_.save(toCString(seqFile)))
        return false;

    return true;
}

inline bool HLAGraphBuilder::load(CharString const & prefix)
{
    // Load GBWT
    CharString gbwtFile = prefix;
    append(gbwtFile, ".gbwt");
    if (!sdsl::load_from_file(gbwt_, toCString(gbwtFile)))
        return false;

    // Load region map
    CharString regionFile = prefix;
    append(regionFile, ".gbwt.regions");
    if (!regionMap_.load(toCString(regionFile)))
        return false;

    // Load coordinate map
    CharString coordFile = prefix;
    append(coordFile, ".gbwt.coords");
    if (!coordMap_.load(toCString(coordFile)))
        return false;

    // Load sequence index
    CharString seqFile = prefix;
    append(seqFile, ".gbwt.seqs");
    if (!seqIndex_.load(toCString(seqFile)))
        return false;

    return true;
}

} // namespace yara_gbwt

#endif // YARA_WITH_GBWT
#endif // APP_YARA_GBWT_HLA_INDEX_H_
