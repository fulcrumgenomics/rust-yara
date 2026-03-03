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
// GSSW-based graph alignment for HLA regions
// ==========================================================================
// This header provides the public interface to GSSW graph alignment.
// The actual GSSW library (which has naming conflicts with yara) is only
// included in the implementation file (gssw_aligner.cpp).

#ifndef APP_YARA_GSSW_ALIGNER_H_
#define APP_YARA_GSSW_ALIGNER_H_

#ifdef YARA_WITH_GBWT

// Include the implementation header which defines the public interface
// using the PIMPL pattern to hide GSSW details
#include "gssw_aligner_impl.h"

// ============================================================================
// Graph building helper functions (for use in indexer)
// ============================================================================

namespace yara_gbwt {

// DNA conversion helper for building graphs
inline char seqanDna5ToChar(uint8_t ord) {
    return dna5ToChar(ord);
}

// Build a GSSW graph from aligned sequences (MSA format)
// Template so it can work with SeqAn StringSets
template <typename TAlleleSeqs>
bool buildGSSWGraphFromMSA(GSSWGraph & graph, TAlleleSeqs const & alleles, uint32_t alleleSetId)
{
    if (empty(alleles)) return false;

    graph.setAlleleSetId(alleleSetId);

    size_t numAlleles = length(alleles);
    size_t msaLength = length(alleles[0]);

    // Verify all sequences have same length (MSA requirement)
    for (size_t i = 1; i < numAlleles; ++i) {
        if (length(alleles[i]) != msaLength) {
            return false;  // Not a valid MSA
        }
    }

    // Find conserved vs variant columns
    std::vector<bool> isConserved(msaLength, true);
    for (size_t col = 0; col < msaLength; ++col) {
        auto first = alleles[0][col];
        for (size_t a = 1; a < numAlleles; ++a) {
            if (alleles[a][col] != first) {
                isConserved[col] = false;
                break;
            }
        }
    }

    // Build nodes from conserved regions and variant positions
    uint32_t nodeId = 0;
    size_t regionStart = 0;
    std::vector<std::vector<uint32_t>> alleleNodes(numAlleles);

    while (regionStart < msaLength) {
        size_t regionEnd = regionStart;

        if (isConserved[regionStart]) {
            while (regionEnd < msaLength && isConserved[regionEnd]) {
                ++regionEnd;
            }

            std::string seq;
            for (size_t i = regionStart; i < regionEnd; ++i) {
                uint8_t ord = ordValue(alleles[0][i]);
                char c = seqanDna5ToChar(ord);
                if (c != '-' && c != 'N') {
                    seq += c;
                }
            }

            if (!seq.empty()) {
                graph.addNode(nodeId, seq, regionStart, regionEnd);
                for (size_t a = 0; a < numAlleles; ++a) {
                    alleleNodes[a].push_back(nodeId);
                }
                ++nodeId;
            }
        } else {
            while (regionEnd < msaLength && !isConserved[regionEnd]) {
                ++regionEnd;
            }

            std::unordered_map<std::string, uint32_t> variantToNode;

            for (size_t a = 0; a < numAlleles; ++a) {
                std::string varSeq;
                for (size_t i = regionStart; i < regionEnd; ++i) {
                    uint8_t ord = ordValue(alleles[a][i]);
                    char c = seqanDna5ToChar(ord);
                    if (c != '-' && c != 'N') {
                        varSeq += c;
                    }
                }

                if (!varSeq.empty()) {
                    auto it = variantToNode.find(varSeq);
                    if (it == variantToNode.end()) {
                        graph.addNode(nodeId, varSeq, regionStart, regionEnd);
                        variantToNode[varSeq] = nodeId;
                        alleleNodes[a].push_back(nodeId);
                        ++nodeId;
                    } else {
                        alleleNodes[a].push_back(it->second);
                    }
                }
            }
        }

        regionStart = regionEnd;
    }

    // Build edges based on allele paths
    for (size_t a = 0; a < numAlleles; ++a) {
        for (size_t i = 1; i < alleleNodes[a].size(); ++i) {
            graph.addEdge(alleleNodes[a][i-1], alleleNodes[a][i]);
        }
    }

    return true;
}

// Build a GSSW graph from unaligned sequences
template <typename TAlleleSeqs>
bool buildGSSWGraphFromSequences(GSSWGraph & graph, TAlleleSeqs const & alleles, uint32_t alleleSetId)
{
    if (empty(alleles)) return false;

    graph.setAlleleSetId(alleleSetId);

    // For unaligned sequences, just use the first (reference) allele
    std::string refSeq;
    for (size_t i = 0; i < length(alleles[0]); ++i) {
        uint8_t ord = ordValue(alleles[0][i]);
        refSeq += seqanDna5ToChar(ord);
    }

    graph.addNode(0, refSeq, 0, length(alleles[0]));

    return true;
}

} // namespace yara_gbwt

#endif // YARA_WITH_GBWT
#endif // APP_YARA_GSSW_ALIGNER_H_
