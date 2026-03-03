// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// GSSW-based graph alignment implementation
// ==========================================================================
// This implementation file includes gssw.h BEFORE any yara headers to avoid
// the name collision between GSSW's 'Match' enum and yara's 'struct Match'.

#ifdef YARA_WITH_GBWT

// Include GSSW first, before any yara headers that define 'Match'
#include <gssw.h>

// Standard includes
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <cstring>
#include <fstream>

// SeqAn includes
#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan2;

namespace yara_gbwt {

// ============================================================================
// DNA conversion helper
// ============================================================================

static const char DNA5_TO_CHAR[] = {'A', 'C', 'G', 'T', 'N'};

// ============================================================================
// GSSWGraphImpl - Internal implementation
// ============================================================================

struct GSSWGraphImpl
{
    struct Node {
        uint32_t id;
        std::string sequence;
        uint32_t refStart;
        uint32_t refEnd;
    };

    struct Edge {
        uint32_t fromId;
        uint32_t toId;
    };

    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::unordered_map<uint32_t, std::vector<uint32_t>> outEdges;
    std::unordered_map<uint32_t, std::vector<uint32_t>> inEdges;
    uint32_t alleleSetId = 0;

    // GSSW structures
    gssw_graph* gsswGraph = nullptr;
    int8_t* nt_table = nullptr;
    int8_t* mat = nullptr;
    std::vector<gssw_node*> gsswNodes;

    ~GSSWGraphImpl() {
        destroy();
    }

    void destroy() {
        if (gsswGraph) {
            gssw_graph_destroy(gsswGraph);
            gsswGraph = nullptr;
        }
        if (nt_table) {
            free(nt_table);
            nt_table = nullptr;
        }
        if (mat) {
            free(mat);
            mat = nullptr;
        }
        gsswNodes.clear();
    }

    void addNode(uint32_t id, std::string const & seq, uint32_t refStart, uint32_t refEnd) {
        nodes.push_back({id, seq, refStart, refEnd});
    }

    void addEdge(uint32_t fromId, uint32_t toId) {
        edges.push_back({fromId, toId});
        outEdges[fromId].push_back(toId);
        inEdges[toId].push_back(fromId);
    }

    bool initGSSW(int8_t match, int8_t mismatch, uint8_t gapOpen, uint8_t gapExtend) {
        destroy();
        if (nodes.empty()) return false;

        nt_table = gssw_create_nt_table();
        if (!nt_table) return false;

        mat = gssw_create_score_matrix(match, mismatch);
        if (!mat) {
            free(nt_table);
            nt_table = nullptr;
            return false;
        }

        gsswGraph = gssw_graph_create(nodes.size());
        if (!gsswGraph) {
            free(nt_table);
            free(mat);
            nt_table = nullptr;
            mat = nullptr;
            return false;
        }

        gsswNodes.resize(nodes.size(), nullptr);
        for (size_t i = 0; i < nodes.size(); ++i) {
            gssw_node* node = gssw_node_create(
                nullptr,
                nodes[i].id,
                nodes[i].sequence.c_str(),
                nt_table,
                mat
            );
            if (!node) {
                destroy();
                return false;
            }
            gsswNodes[i] = node;
            gssw_graph_add_node(gsswGraph, node);
        }

        for (auto const & edge : edges) {
            if (edge.fromId < gsswNodes.size() && edge.toId < gsswNodes.size()) {
                gssw_nodes_add_edge(gsswNodes[edge.fromId], gsswNodes[edge.toId]);
            }
        }

        return true;
    }

    int32_t align(const char* readSeq, int32_t readLen, uint32_t & matchBegin,
                  uint32_t & matchEnd, uint32_t & errors) {
        if (!gsswGraph || !nt_table || !mat) {
            return -1;
        }

        gssw_graph_fill(gsswGraph, readSeq, nt_table, mat,
                        3,         // gap open
                        1,         // gap extend
                        0,         // start_full_length_bonus
                        0,         // end_full_length_bonus
                        readLen,   // maskLen
                        2,         // score_size
                        false      // save_matrixes
        );

        gssw_graph_mapping* mapping = gssw_graph_trace_back(
            gsswGraph,
            readSeq,
            readLen,
            nt_table,
            mat,
            3,   // gap open
            1,   // gap extend
            0,   // start_full_length_bonus
            0    // end_full_length_bonus
        );

        if (!mapping) {
            return -1;
        }

        int32_t score = mapping->score;

        errors = 0;
        matchBegin = 0;
        matchEnd = 0;

        if (mapping->cigar.length > 0) {
            gssw_node_cigar* firstNC = &mapping->cigar.elements[0];
            gssw_node_cigar* lastNC = &mapping->cigar.elements[mapping->cigar.length - 1];

            for (size_t i = 0; i < nodes.size(); ++i) {
                if (nodes[i].id == firstNC->node->id) {
                    matchBegin = nodes[i].refStart;
                }
                if (nodes[i].id == lastNC->node->id) {
                    matchEnd = nodes[i].refEnd;
                }
            }

            for (int32_t i = 0; i < mapping->cigar.length; ++i) {
                gssw_node_cigar* nc = &mapping->cigar.elements[i];
                for (int32_t j = 0; j < nc->cigar->length; ++j) {
                    char op = nc->cigar->elements[j].type;
                    uint32_t len = nc->cigar->elements[j].length;
                    if (op == 'X' || op == 'I' || op == 'D') {
                        errors += len;
                    }
                }
            }
        }

        gssw_graph_mapping_destroy(mapping);
        return score;
    }
};

} // namespace yara_gbwt

// Now include the header which defines the public interface
#include "gssw_aligner_impl.h"

namespace yara_gbwt {

// ============================================================================
// GSSWGraph implementation (delegating to impl)
// ============================================================================

GSSWGraph::GSSWGraph() : impl_(std::make_unique<GSSWGraphImpl>()) {}
GSSWGraph::~GSSWGraph() = default;

void GSSWGraph::addNode(uint32_t id, std::string const & seq, uint32_t refStart, uint32_t refEnd) {
    impl_->addNode(id, seq, refStart, refEnd);
}

void GSSWGraph::addEdge(uint32_t fromId, uint32_t toId) {
    impl_->addEdge(fromId, toId);
}

bool GSSWGraph::initGSSW(int8_t match, int8_t mismatch, uint8_t gapOpen, uint8_t gapExtend) {
    return impl_->initGSSW(match, mismatch, gapOpen, gapExtend);
}

size_t GSSWGraph::nodeCount() const { return impl_->nodes.size(); }
size_t GSSWGraph::edgeCount() const { return impl_->edges.size(); }
uint32_t GSSWGraph::alleleSetId() const { return impl_->alleleSetId; }
void GSSWGraph::setAlleleSetId(uint32_t id) { impl_->alleleSetId = id; }

int32_t GSSWGraph::alignImpl(const char* readSeq, int32_t readLen, uint32_t & matchBegin,
                              uint32_t & matchEnd, uint32_t & errors) {
    return impl_->align(readSeq, readLen, matchBegin, matchEnd, errors);
}

bool GSSWGraph::getRefPosition(uint32_t nodeId, uint32_t posInNode, uint32_t & refPos) const {
    for (auto const & node : impl_->nodes) {
        if (node.id == nodeId) {
            refPos = node.refStart + posInNode;
            return true;
        }
    }
    return false;
}

bool GSSWGraph::save(const char * fileName) const {
    std::ofstream out(fileName, std::ios::binary);
    if (!out) return false;

    out.write(reinterpret_cast<const char*>(&impl_->alleleSetId), sizeof(impl_->alleleSetId));

    size_t numNodes = impl_->nodes.size();
    out.write(reinterpret_cast<const char*>(&numNodes), sizeof(numNodes));
    for (auto const & node : impl_->nodes) {
        out.write(reinterpret_cast<const char*>(&node.id), sizeof(node.id));
        size_t seqLen = node.sequence.size();
        out.write(reinterpret_cast<const char*>(&seqLen), sizeof(seqLen));
        out.write(node.sequence.data(), seqLen);
        out.write(reinterpret_cast<const char*>(&node.refStart), sizeof(node.refStart));
        out.write(reinterpret_cast<const char*>(&node.refEnd), sizeof(node.refEnd));
    }

    size_t numEdges = impl_->edges.size();
    out.write(reinterpret_cast<const char*>(&numEdges), sizeof(numEdges));
    for (auto const & edge : impl_->edges) {
        out.write(reinterpret_cast<const char*>(&edge.fromId), sizeof(edge.fromId));
        out.write(reinterpret_cast<const char*>(&edge.toId), sizeof(edge.toId));
    }

    return out.good();
}

bool GSSWGraph::load(const char * fileName) {
    std::ifstream in(fileName, std::ios::binary);
    if (!in) return false;

    impl_->destroy();
    impl_->nodes.clear();
    impl_->edges.clear();
    impl_->outEdges.clear();
    impl_->inEdges.clear();

    in.read(reinterpret_cast<char*>(&impl_->alleleSetId), sizeof(impl_->alleleSetId));

    size_t numNodes;
    in.read(reinterpret_cast<char*>(&numNodes), sizeof(numNodes));
    for (size_t i = 0; i < numNodes; ++i) {
        GSSWGraphImpl::Node node;
        in.read(reinterpret_cast<char*>(&node.id), sizeof(node.id));
        size_t seqLen;
        in.read(reinterpret_cast<char*>(&seqLen), sizeof(seqLen));
        node.sequence.resize(seqLen);
        in.read(&node.sequence[0], seqLen);
        in.read(reinterpret_cast<char*>(&node.refStart), sizeof(node.refStart));
        in.read(reinterpret_cast<char*>(&node.refEnd), sizeof(node.refEnd));
        impl_->nodes.push_back(node);
    }

    size_t numEdges;
    in.read(reinterpret_cast<char*>(&numEdges), sizeof(numEdges));
    for (size_t i = 0; i < numEdges; ++i) {
        GSSWGraphImpl::Edge edge;
        in.read(reinterpret_cast<char*>(&edge.fromId), sizeof(edge.fromId));
        in.read(reinterpret_cast<char*>(&edge.toId), sizeof(edge.toId));
        impl_->edges.push_back(edge);
        impl_->outEdges[edge.fromId].push_back(edge.toId);
        impl_->inEdges[edge.toId].push_back(edge.fromId);
    }

    return in.good();
}

// ============================================================================
// GSSWAligner implementation
// ============================================================================

bool GSSWAligner::save(const char * prefix) const {
    std::string indexFile = std::string(prefix) + ".gssw.index";
    std::ofstream out(indexFile, std::ios::binary);
    if (!out) return false;

    size_t numGraphs = graphs_.size();
    out.write(reinterpret_cast<const char*>(&numGraphs), sizeof(numGraphs));

    for (auto const & pair : graphs_) {
        uint32_t setId = pair.first;
        out.write(reinterpret_cast<const char*>(&setId), sizeof(setId));

        std::string graphFile = std::string(prefix) + ".gssw." + std::to_string(setId);
        if (!pair.second->save(graphFile.c_str())) {
            return false;
        }
    }

    return out.good();
}

bool GSSWAligner::load(const char * prefix) {
    graphs_.clear();

    std::string indexFile = std::string(prefix) + ".gssw.index";
    std::ifstream in(indexFile, std::ios::binary);
    if (!in) return false;

    size_t numGraphs;
    in.read(reinterpret_cast<char*>(&numGraphs), sizeof(numGraphs));

    for (size_t i = 0; i < numGraphs; ++i) {
        uint32_t setId;
        in.read(reinterpret_cast<char*>(&setId), sizeof(setId));

        auto graph = std::make_shared<GSSWGraph>();
        std::string graphFile = std::string(prefix) + ".gssw." + std::to_string(setId);
        if (!graph->load(graphFile.c_str())) {
            return false;
        }
        if (!graph->initGSSW()) {
            return false;
        }

        graphs_[setId] = graph;
    }

    return true;
}

} // namespace yara_gbwt

#endif // YARA_WITH_GBWT
