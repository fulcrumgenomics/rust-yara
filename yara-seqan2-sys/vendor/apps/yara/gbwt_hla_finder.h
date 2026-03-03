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
// GBWT-based seed finding for HLA regions
// ==========================================================================

#ifndef APP_YARA_GBWT_HLA_FINDER_H_
#define APP_YARA_GBWT_HLA_FINDER_H_

#ifdef YARA_WITH_GBWT

#include <gbwt/gbwt.h>
#include "gbwt_hla_index.h"

using namespace seqan2;

namespace yara_gbwt {

// ============================================================================
// GBWTSeedFinder - Find seeds in HLA regions using sequence index
// ============================================================================
// Uses SeqAn's FM-index over concatenated allele sequences to find seed
// occurrences across all HLA alleles simultaneously.

template <typename TReadSeq>
class GBWTSeedFinder
{
public:
    typedef typename Size<TReadSeq>::Type TSize;

    // Result structure for a found seed
    struct SeedHit {
        uint32_t alleleId;          // Which allele was matched
        uint32_t alleleSetId;       // Which HLA gene (set of alleles)
        uint32_t posInAllele;       // Position within the allele sequence
        TSize seedErrors;           // Errors in seed (0 for exact)
        TSize seedBegin;            // Seed start in read
        TSize seedEnd;              // Seed end in read
    };

    GBWTSeedFinder(HLASequenceIndex const & seqIndex,
                   AlleleCoordMap const & coordMap,
                   HLARegionMap const & regionMap) :
        seqIndex_(seqIndex), coordMap_(coordMap), regionMap_(regionMap)
    {}

    // Find exact seed matches in HLA allele sequences
    // Reports all matching positions via delegate
    template <typename TDelegate>
    void findSeeds(TReadSeq const & read,
                   TSize seedLength,
                   TSize seedStep,
                   TDelegate && delegate);

    // Find seeds with up to maxErrors mismatches
    template <typename TDelegate>
    void findSeedsApprox(TReadSeq const & read,
                         TSize seedLength,
                         TSize seedStep,
                         TSize maxErrors,
                         TDelegate && delegate);

private:
    HLASequenceIndex const & seqIndex_;
    AlleleCoordMap const & coordMap_;
    HLARegionMap const & regionMap_;

    // Process FM-index occurrences and convert to SeedHits
    template <typename TDelegate, typename TOccurrences>
    void processOccurrences_(TOccurrences const & occs,
                             TSize seedBegin,
                             TSize seedEnd,
                             TSize seedErrors,
                             TDelegate && delegate);
};

// ============================================================================
// Implementation
// ============================================================================

template <typename TReadSeq>
template <typename TDelegate>
void GBWTSeedFinder<TReadSeq>::findSeeds(
    TReadSeq const & read,
    TSize seedLength,
    TSize seedStep,
    TDelegate && delegate)
{
    TSize readLen = length(read);
    if (readLen < seedLength) return;

    typedef typename HLASequenceIndex::TIndex TIndex;
    typedef Finder<TIndex> TFinder;

    // Extract seeds at regular intervals
    for (TSize seedStart = 0; seedStart + seedLength <= readLen; seedStart += seedStep)
    {
        // Get seed sequence
        auto seed = infix(read, seedStart, seedStart + seedLength);

        // Use FM-index exact search
        TFinder finder(seqIndex_.index);

        // Find exact occurrences
        while (find(finder, seed))
        {
            // Get the position in the concatenated sequence
            auto pos = position(finder);

            // pos is a Pair<seqId, posInSeq>
            uint32_t alleleId = getSeqNo(pos);
            uint32_t posInAllele = getSeqOffset(pos);

            if (alleleId < seqIndex_.alleleSetIds.size())
            {
                SeedHit hit;
                hit.alleleId = alleleId;
                hit.alleleSetId = seqIndex_.alleleSetIds[alleleId];
                hit.posInAllele = posInAllele;
                hit.seedErrors = 0;
                hit.seedBegin = seedStart;
                hit.seedEnd = seedStart + seedLength;

                delegate(hit);
            }
        }
    }
}

template <typename TReadSeq>
template <typename TDelegate>
void GBWTSeedFinder<TReadSeq>::findSeedsApprox(
    TReadSeq const & read,
    TSize seedLength,
    TSize seedStep,
    TSize maxErrors,
    TDelegate && delegate)
{
    TSize readLen = length(read);
    if (readLen < seedLength) return;

    typedef typename HLASequenceIndex::TIndex TIndex;

    for (TSize seedStart = 0; seedStart + seedLength <= readLen; seedStart += seedStep)
    {
        auto seed = infix(read, seedStart, seedStart + seedLength);

        // For approximate search, we iterate through alleles and check distance
        // This is a simple implementation; could be optimized with backtracking search
        for (size_t alleleId = 0; alleleId < length(seqIndex_.alleleSeqs); ++alleleId)
        {
            auto const & alleleSeq = seqIndex_.alleleSeqs[alleleId];
            if (length(alleleSeq) < seedLength) continue;

            // Slide through the allele looking for approximate matches
            for (TSize pos = 0; pos + seedLength <= length(alleleSeq); ++pos)
            {
                auto alleleWindow = infix(alleleSeq, pos, pos + seedLength);

                // Count mismatches
                TSize errors = 0;
                for (TSize i = 0; i < seedLength && errors <= maxErrors; ++i)
                {
                    if (seed[i] != alleleWindow[i])
                        ++errors;
                }

                if (errors <= maxErrors)
                {
                    SeedHit hit;
                    hit.alleleId = alleleId;
                    hit.alleleSetId = seqIndex_.alleleSetIds[alleleId];
                    hit.posInAllele = pos;
                    hit.seedErrors = errors;
                    hit.seedBegin = seedStart;
                    hit.seedEnd = seedStart + seedLength;

                    delegate(hit);
                }
            }
        }
    }
}

} // namespace yara_gbwt

#endif // YARA_WITH_GBWT
#endif // APP_YARA_GBWT_HLA_FINDER_H_
