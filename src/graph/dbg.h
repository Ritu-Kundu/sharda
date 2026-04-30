#pragma once

#include "graph/types.h"
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

namespace sharda {

/// Positional De Bruijn Graph.
/// Backbone nodes are indexed by ref_pos; non-backbone nodes by kmer string.
class DBG {
public:
    explicit DBG(int k) : k_(k) {}

    int k() const { return k_; }

    // ── Node operations ─────────────────────────────────────────────────
    /// Add a backbone node at the given reference position.
    uint64_t add_backbone_node(const std::string& kmer, int32_t ref_pos, int tr_id);

    /// Find or create a non-backbone node for a kmer.
    uint64_t add_read_node(const std::string& kmer);

    /// Get node by ID.
    const Node& node(uint64_t id) const { return nodes_[id]; }
    Node& node_mut(uint64_t id) { return nodes_[id]; }

    /// All nodes.
    const std::vector<Node>& nodes() const { return nodes_; }
    size_t node_count() const { return nodes_.size(); }
    size_t active_node_count() const;
    bool is_node_removed(uint64_t node_id) const;

    /// Lookup backbone node by ref_pos.  Returns UINT64_MAX if not found.
    uint64_t backbone_node_at(int32_t ref_pos) const;

    /// Lookup non-backbone node by kmer. Returns UINT64_MAX if not found.
    uint64_t find_read_node(const std::string& kmer) const;

    // ── Edge operations ─────────────────────────────────────────────────
    /// Add or increment a directed edge.
    void add_edge(uint64_t from, uint64_t to);

    /// All edges (as flat list).
    const std::vector<Edge>& edges() const { return edges_; }
    size_t edge_count() const { return edges_.size(); }

    /// Get outgoing edges for a node.
    const std::vector<uint64_t>& out_edges(uint64_t node_id) const;

    /// Get incoming edges for a node.
    const std::vector<uint64_t>& in_edges(uint64_t node_id) const;

    // ── Haplotype edges ─────────────────────────────────────────────────
    void add_haplotype_edge(uint64_t from, uint64_t to);
    const std::vector<HaplotypeEdge>& haplotype_edges() const { return hap_edges_; }

    // ── Backbone TR helpers ─────────────────────────────────────────────
    const std::vector<uint64_t>& tr_nodes(int tr_id) const;

    // ── Mutation for cleaning ───────────────────────────────────────────
    void remove_node(uint64_t id);
    void remove_edge(uint64_t edge_idx);

    /// Rebuild adjacency lists after removals.
    void rebuild_adjacency();

private:
    int k_;

    std::vector<Node>  nodes_;
    std::vector<Edge>  edges_;
    std::vector<HaplotypeEdge> hap_edges_;

    // Adjacency: node_id → list of edge indices
    std::vector<std::vector<uint64_t>> out_adj_;
    std::vector<std::vector<uint64_t>> in_adj_;

    // Lookup maps
    std::unordered_map<int32_t, uint64_t>     pos_to_node_;   // ref_pos → node_id
    std::unordered_map<std::string, uint64_t>  kmer_to_node_;  // kmer → node_id (non-backbone)
    std::unordered_map<int, std::vector<uint64_t>> tr_to_nodes_; // tr_id → backbone node ids

    // edge lookup: (from,to) → edge index
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, uint64_t>> edge_map_;

    // Removed flags
    std::vector<bool> node_removed_;
    std::vector<bool> edge_removed_;

    static const std::vector<uint64_t> empty_vec_;

    void ensure_adj_size(uint64_t id);
};

} // namespace sharda
