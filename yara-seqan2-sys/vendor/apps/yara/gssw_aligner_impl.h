// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// GSSW aligner implementation header - internal use only
// ==========================================================================

#ifndef APP_YARA_GSSW_ALIGNER_IMPL_H_
#define APP_YARA_GSSW_ALIGNER_IMPL_H_

#ifdef YARA_WITH_GBWT

#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <cstdint>

namespace yara_gbwt {

// Forward declaration
struct GSSWGraphImpl;

// ============================================================================
// GSSWGraph - Public interface (PIMPL pattern)
// ============================================================================

class GSSWGraph
{
public:
    GSSWGraph();
    ~GSSWGraph();

    // Disable copy, allow move
    GSSWGraph(const GSSWGraph&) = delete;
    GSSWGraph& operator=(const GSSWGraph&) = delete;
    GSSWGraph(GSSWGraph&&) = default;
    GSSWGraph& operator=(GSSWGraph&&) = default;

    void addNode(uint32_t id, std::string const & seq, uint32_t refStart, uint32_t refEnd);
    void addEdge(uint32_t fromId, uint32_t toId);
    bool initGSSW(int8_t match = 2, int8_t mismatch = -2,
                  uint8_t gapOpen = 3, uint8_t gapExtend = 1);

    // Template method that converts read to string and calls impl
    template <typename TReadSeq>
    int32_t align(TReadSeq const & read, uint32_t & matchBegin, uint32_t & matchEnd,
                  uint32_t & errors);

    bool getRefPosition(uint32_t nodeId, uint32_t posInNode, uint32_t & refPos) const;

    bool save(const char * fileName) const;
    bool load(const char * fileName);

    size_t nodeCount() const;
    size_t edgeCount() const;
    uint32_t alleleSetId() const;
    void setAlleleSetId(uint32_t id);

    // Internal: called by template method
    int32_t alignImpl(const char* readSeq, int32_t readLen, uint32_t & matchBegin,
                      uint32_t & matchEnd, uint32_t & errors);

private:
    std::unique_ptr<GSSWGraphImpl> impl_;
};

// ============================================================================
// GSSWAligner - Manages multiple graphs
// ============================================================================

class GSSWAligner
{
public:
    struct AlignmentResult {
        int32_t score;
        uint32_t refBegin;
        uint32_t refEnd;
        uint32_t errors;
        bool valid;
    };

    GSSWAligner() {}

    void addGraph(uint32_t alleleSetId, std::shared_ptr<GSSWGraph> graph) {
        graphs_[alleleSetId] = graph;
    }

    template <typename TReadSeq>
    AlignmentResult align(TReadSeq const & read, uint32_t alleleSetId);

    bool hasGraph(uint32_t alleleSetId) const {
        return graphs_.find(alleleSetId) != graphs_.end();
    }

    bool save(const char * prefix) const;
    bool load(const char * prefix);

private:
    std::unordered_map<uint32_t, std::shared_ptr<GSSWGraph>> graphs_;
};

// ============================================================================
// Template implementations
// ============================================================================

// DNA conversion - must be inline since it's in template code
inline char dna5ToChar(uint8_t ord) {
    static const char DNA5_TO_CHAR[] = {'A', 'C', 'G', 'T', 'N'};
    return (ord < 5) ? DNA5_TO_CHAR[ord] : 'N';
}

template <typename TReadSeq>
int32_t GSSWGraph::align(TReadSeq const & read, uint32_t & matchBegin,
                          uint32_t & matchEnd, uint32_t & errors)
{
    std::string readStr;
    readStr.reserve(length(read));
    for (size_t i = 0; i < length(read); ++i) {
        readStr += dna5ToChar(ordValue(read[i]));
    }

    return alignImpl(readStr.c_str(), static_cast<int32_t>(readStr.length()),
                     matchBegin, matchEnd, errors);
}

template <typename TReadSeq>
GSSWAligner::AlignmentResult GSSWAligner::align(TReadSeq const & read, uint32_t alleleSetId)
{
    AlignmentResult result = {0, 0, 0, 0, false};

    auto it = graphs_.find(alleleSetId);
    if (it == graphs_.end() || !it->second) {
        return result;
    }

    int32_t score = it->second->align(read, result.refBegin, result.refEnd, result.errors);
    if (score > 0) {
        result.score = score;
        result.valid = true;
    }

    return result;
}

} // namespace yara_gbwt

#endif // YARA_WITH_GBWT
#endif // APP_YARA_GSSW_ALIGNER_IMPL_H_
