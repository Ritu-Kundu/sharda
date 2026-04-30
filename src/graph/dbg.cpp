#include "graph/dbg.h"
#include <spdlog/spdlog.h>
#include <algorithm>

namespace sharda {

const std::vector<uint64_t> DBG::empty_vec_ = {};

void DBG::ensure_adj_size(uint64_t id) {
    if (id >= out_adj_.size()) {
        out_adj_.resize(id + 1);
        in_adj_.resize(id + 1);
        node_removed_.resize(id + 1, false);
    }
}

uint64_t DBG::add_backbone_node(const std::string& kmer, int32_t ref_pos, int tr_id) {
    uint64_t id = nodes_.size();
    nodes_.push_back({id, kmer, ref_pos, tr_id, true, 0});
    pos_to_node_[ref_pos] = id;
    if (tr_id >= 0) {
        tr_to_nodes_[tr_id].push_back(id);
    }
    ensure_adj_size(id);
    return id;
}

uint64_t DBG::add_read_node(const std::string& kmer) {
    auto it = kmer_to_node_.find(kmer);
    if (it != kmer_to_node_.end()) return it->second;

    uint64_t id = nodes_.size();
    nodes_.push_back({id, kmer, -1, -1, false, 0});
    kmer_to_node_[kmer] = id;
    ensure_adj_size(id);
    return id;
}

uint64_t DBG::backbone_node_at(int32_t ref_pos) const {
    auto it = pos_to_node_.find(ref_pos);
    return (it != pos_to_node_.end()) ? it->second : UINT64_MAX;
}

uint64_t DBG::find_read_node(const std::string& kmer) const {
    auto it = kmer_to_node_.find(kmer);
    return (it != kmer_to_node_.end()) ? it->second : UINT64_MAX;
}

size_t DBG::active_node_count() const {
    size_t active = 0;
    for (size_t i = 0; i < nodes_.size(); ++i) {
        if (!is_node_removed(i)) {
            active++;
        }
    }
    return active;
}

bool DBG::is_node_removed(uint64_t node_id) const {
    return node_id < node_removed_.size() && node_removed_[node_id];
}

void DBG::add_edge(uint64_t from, uint64_t to) {
    auto it = edge_map_.find(from);
    if (it != edge_map_.end()) {
        auto jt = it->second.find(to);
        if (jt != it->second.end()) {
            edges_[jt->second].weight++;
            return;
        }
    }
    uint64_t idx = edges_.size();
    edges_.push_back({from, to, 1});
    edge_removed_.push_back(false);
    edge_map_[from][to] = idx;

    ensure_adj_size(std::max(from, to));
    out_adj_[from].push_back(idx);
    in_adj_[to].push_back(idx);
}

const std::vector<uint64_t>& DBG::out_edges(uint64_t node_id) const {
    return (node_id < out_adj_.size()) ? out_adj_[node_id] : empty_vec_;
}

const std::vector<uint64_t>& DBG::in_edges(uint64_t node_id) const {
    return (node_id < in_adj_.size()) ? in_adj_[node_id] : empty_vec_;
}

void DBG::add_haplotype_edge(uint64_t from, uint64_t to) {
    hap_edges_.push_back({from, to, 1});
}

const std::vector<uint64_t>& DBG::tr_nodes(int tr_id) const {
    auto it = tr_to_nodes_.find(tr_id);
    return (it != tr_to_nodes_.end()) ? it->second : empty_vec_;
}

void DBG::remove_node(uint64_t id) {
    if (id < node_removed_.size()) node_removed_[id] = true;
}

void DBG::remove_edge(uint64_t edge_idx) {
    if (edge_idx < edge_removed_.size()) edge_removed_[edge_idx] = true;
}

void DBG::rebuild_adjacency() {
    // Rebuild out_adj_ and in_adj_ from non-removed edges
    for (auto& v : out_adj_) v.clear();
    for (auto& v : in_adj_)  v.clear();
    edge_map_.clear();

    std::vector<Edge> new_edges;
    new_edges.reserve(edges_.size());

    for (size_t i = 0; i < edges_.size(); ++i) {
        if (edge_removed_[i]) continue;
        const auto& e = edges_[i];
        if (node_removed_[e.from] || node_removed_[e.to]) continue;
        uint64_t idx = new_edges.size();
        new_edges.push_back(e);
        edge_map_[e.from][e.to] = idx;
        out_adj_[e.from].push_back(idx);
        in_adj_[e.to].push_back(idx);
    }
    edges_ = std::move(new_edges);
    edge_removed_.assign(edges_.size(), false);

    // Also prune haplotype edges referencing removed nodes
    std::vector<HaplotypeEdge> new_hap;
    for (const auto& he : hap_edges_) {
        if (!node_removed_[he.from_node] && !node_removed_[he.to_node])
            new_hap.push_back(he);
    }
    hap_edges_ = std::move(new_hap);

    spdlog::debug("rebuild_adjacency: {} nodes active, {} edges, {} hap_edges",
                  nodes_.size(), edges_.size(), hap_edges_.size());
}

} // namespace sharda
