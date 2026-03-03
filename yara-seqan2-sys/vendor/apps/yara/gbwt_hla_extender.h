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
// GSSW-based graph extension for HLA regions
// ==========================================================================
// Uses GSSW (Graph Smith-Waterman) to align reads to variation graphs
// built from aligned HLA allele sequences. This provides O(read_length × graph_width)
// complexity instead of O(read_length × num_alleles).

#ifndef APP_YARA_GBWT_HLA_EXTENDER_H_
#define APP_YARA_GBWT_HLA_EXTENDER_H_

#ifdef YARA_WITH_GBWT

#include <gbwt/gbwt.h>
#include "gbwt_hla_index.h"
#include "gbwt_hla_finder.h"
#include "gssw_aligner.h"
#include "find_extender.h"

using namespace seqan2;

namespace yara_gbwt {

// ============================================================================
// GBWTHLAExtender - Extend seeds through HLA variation graphs using GSSW
// ============================================================================
// Uses GSSW graph alignment for HLA regions when a graph is available,
// falling back to per-allele Myers extension otherwise.

template <typename TContigSeqs, typename TReadSeq>
class GBWTHLAExtender
{
public:
    typedef typename Size<TReadSeq>::Type TSize;
    typedef typename StringSetPosition<TContigSeqs>::Type TContigsPos;

    // Allele sequence types for fallback
    typedef typename HLASequenceIndex::TAlleleSeqs TAlleleSeqs;
    typedef typename StringSetPosition<TAlleleSeqs>::Type TAllelePos;
    typedef Extender<TAlleleSeqs, TReadSeq> TAlleleExtender;

    struct ExtensionResult {
        TContigsPos refBegin;       // Start position in reference
        TContigsPos refEnd;         // End position in reference
        TSize matchBegin;           // Start in read
        TSize matchEnd;             // End in read
        TSize errors;               // Total errors
        uint32_t alleleId;          // Which allele this match is from (0 if from graph)
        bool valid;                 // Whether extension succeeded
    };

    GBWTHLAExtender(HLASequenceIndex const & seqIndex,
                    AlleleCoordMap const & coordMap,
                    HLARegionMap const & regionMap,
                    GSSWAligner * gsswAligner = nullptr) :
        seqIndex_(seqIndex),
        coordMap_(coordMap),
        regionMap_(regionMap),
        gsswAligner_(gsswAligner),
        alleleExtender_(seqIndex.alleleSeqs)
    {}

    // Extend from a seed hit using GSSW graph alignment
    template <typename TDistance, typename TErrors, typename TDelegate>
    void extendFromSeed(TReadSeq const & read,
                        typename GBWTSeedFinder<TReadSeq>::SeedHit const & seedHit,
                        TDistance const & /* tag */,
                        TErrors maxErrors,
                        TDelegate && delegate);

    // Align entire read to graph (for rescue path)
    template <typename TErrors, typename TDelegate>
    void alignToGraph(TReadSeq const & read,
                      uint32_t alleleSetId,
                      uint32_t refContigId,
                      TErrors maxErrors,
                      TDelegate && delegate);

    // Access for external use
    HLASequenceIndex const & seqIndex() const { return seqIndex_; }
    AlleleCoordMap const & coordMap() const { return coordMap_; }
    HLARegionMap const & regionMap() const { return regionMap_; }

private:
    HLASequenceIndex const & seqIndex_;
    AlleleCoordMap const & coordMap_;
    HLARegionMap const & regionMap_;
    GSSWAligner * gsswAligner_;
    TAlleleExtender alleleExtender_;

    // Fallback to per-allele Myers extension
    template <typename TDistance, typename TErrors, typename TDelegate>
    void extendWithMyers_(TReadSeq const & read,
                          typename GBWTSeedFinder<TReadSeq>::SeedHit const & seedHit,
                          TDistance const &,
                          TErrors maxErrors,
                          TDelegate && delegate);
};

// ============================================================================
// Implementation
// ============================================================================

template <typename TContigSeqs, typename TReadSeq>
template <typename TDistance, typename TErrors, typename TDelegate>
void GBWTHLAExtender<TContigSeqs, TReadSeq>::extendFromSeed(
    TReadSeq const & read,
    typename GBWTSeedFinder<TReadSeq>::SeedHit const & seedHit,
    TDistance const & distTag,
    TErrors maxErrors,
    TDelegate && delegate)
{
    TSize readLen = length(read);
    if (readLen == 0) return;

    // Try GSSW graph alignment first if available
    if (gsswAligner_ && gsswAligner_->hasGraph(seedHit.alleleSetId))
    {
        auto result = gsswAligner_->align(read, seedHit.alleleSetId);
        if (result.valid && result.errors <= static_cast<uint32_t>(maxErrors))
        {
            // Get reference contig ID from coordMap
            auto const * alleleIndices = coordMap_.getAllelesForSet(seedHit.alleleSetId);
            if (alleleIndices && !alleleIndices->empty())
            {
                // Use first allele's contig ID (all alleles in set map to same contig)
                auto const & alleleInfo = coordMap_.alleles[(*alleleIndices)[0]];

                ExtensionResult extResult;
                extResult.valid = true;
                extResult.matchBegin = 0;
                extResult.matchEnd = readLen;
                extResult.errors = result.errors;
                extResult.alleleId = 0;  // Graph-based, not specific allele

                // Map graph coordinates to reference coordinates
                setValueI1(extResult.refBegin, alleleInfo.refContigId);
                setSeqOffset(extResult.refBegin, alleleInfo.refStart + result.refBegin);
                setValueI1(extResult.refEnd, alleleInfo.refContigId);
                setSeqOffset(extResult.refEnd, alleleInfo.refStart + result.refEnd);

                delegate(extResult);
                return;
            }
        }
    }

    // Fall back to per-allele Myers extension
    extendWithMyers_(read, seedHit, distTag, maxErrors, std::forward<TDelegate>(delegate));
}

template <typename TContigSeqs, typename TReadSeq>
template <typename TErrors, typename TDelegate>
void GBWTHLAExtender<TContigSeqs, TReadSeq>::alignToGraph(
    TReadSeq const & read,
    uint32_t alleleSetId,
    uint32_t refContigId,
    TErrors maxErrors,
    TDelegate && delegate)
{
    TSize readLen = length(read);
    if (readLen == 0) return;

    // Try GSSW graph alignment
    if (gsswAligner_ && gsswAligner_->hasGraph(alleleSetId))
    {
        auto result = gsswAligner_->align(read, alleleSetId);
        if (result.valid && result.errors <= static_cast<uint32_t>(maxErrors))
        {
            auto const * alleleIndices = coordMap_.getAllelesForSet(alleleSetId);
            if (alleleIndices && !alleleIndices->empty())
            {
                auto const & alleleInfo = coordMap_.alleles[(*alleleIndices)[0]];

                ExtensionResult extResult;
                extResult.valid = true;
                extResult.matchBegin = 0;
                extResult.matchEnd = readLen;
                extResult.errors = result.errors;
                extResult.alleleId = 0;

                setValueI1(extResult.refBegin, refContigId);
                setSeqOffset(extResult.refBegin, alleleInfo.refStart + result.refBegin);
                setValueI1(extResult.refEnd, refContigId);
                setSeqOffset(extResult.refEnd, alleleInfo.refStart + result.refEnd);

                delegate(extResult);
            }
        }
    }
}

template <typename TContigSeqs, typename TReadSeq>
template <typename TDistance, typename TErrors, typename TDelegate>
void GBWTHLAExtender<TContigSeqs, TReadSeq>::extendWithMyers_(
    TReadSeq const & read,
    typename GBWTSeedFinder<TReadSeq>::SeedHit const & seedHit,
    TDistance const &,
    TErrors maxErrors,
    TDelegate && delegate)
{
    TSize readLen = length(read);

    // Validate allele ID
    if (seedHit.alleleId >= length(seqIndex_.alleleSeqs))
        return;

    TSize alleleLen = length(seqIndex_.alleleSeqs[seedHit.alleleId]);
    if (alleleLen == 0) return;
    if (seedHit.posInAllele >= alleleLen) return;

    TSize seedLen = seedHit.seedEnd - seedHit.seedBegin;
    if (seedLen == 0 || seedLen > readLen) return;

    // Ensure seed doesn't extend past allele end
    TSize seedEndInAllele = seedHit.posInAllele + seedLen;
    if (seedEndInAllele > alleleLen)
        seedEndInAllele = alleleLen;

    // Create allele positions for extension
    TAllelePos alleleBegin;
    setValueI1(alleleBegin, seedHit.alleleId);
    setSeqOffset(alleleBegin, seedHit.posInAllele);

    TAllelePos alleleEnd;
    setValueI1(alleleEnd, seedHit.alleleId);
    setSeqOffset(alleleEnd, seedEndInAllele);

    uint32_t alleleId = seedHit.alleleId;
    AlleleCoordMap const & coordMap = coordMap_;

    // Use Myers extension on allele sequence
    extend(alleleExtender_,
           read,
           alleleBegin, alleleEnd,
           seedHit.seedBegin, seedHit.seedEnd,
           TDistance(),
           seedHit.seedErrors, maxErrors,
           [&, alleleId, readLen](TAllelePos matchBeginAllele, TAllelePos matchEndAllele, TErrors errors) {
               // Verify this allele matches the expected ID
               if (getSeqNo(matchBeginAllele) != alleleId || getSeqNo(matchEndAllele) != alleleId)
                   return;

               // Translate allele coordinates to reference coordinates
               uint32_t allelePosBegin = getSeqOffset(matchBeginAllele);
               uint32_t allelePosEnd = getSeqOffset(matchEndAllele);

               uint32_t refContigId, refPosBegin;
               bool validBegin = coordMap.toRefCoord(alleleId, allelePosBegin, refContigId, refPosBegin);

               uint32_t refContigIdEnd, refPosEnd;
               bool validEnd = coordMap.toRefCoord(alleleId, allelePosEnd > 0 ? allelePosEnd - 1 : 0,
                                                   refContigIdEnd, refPosEnd);

               // Ensure both positions are on the same contig
               if (validBegin && validEnd && refContigId == refContigIdEnd)
               {
                   ExtensionResult result;
                   result.valid = true;
                   result.matchBegin = 0;
                   result.matchEnd = readLen;
                   result.errors = errors;
                   result.alleleId = alleleId;

                   setValueI1(result.refBegin, refContigId);
                   setSeqOffset(result.refBegin, refPosBegin);
                   setValueI1(result.refEnd, refContigId);
                   setSeqOffset(result.refEnd, refPosEnd + 1);

                   delegate(result);
               }
           });
}

} // namespace yara_gbwt

#endif // YARA_WITH_GBWT
#endif // APP_YARA_GBWT_HLA_EXTENDER_H_
