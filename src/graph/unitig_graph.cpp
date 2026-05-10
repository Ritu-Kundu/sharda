#include "graph/unitig_graph.h"
#include <spdlog/spdlog.h>
#include <algorithm>
#include <unordered_set>
#include <stack>
#include <cmath>

namespace sharda {

uint64_t UnitigGraph::node_to_unitig(uint64_t node_id) const {
    auto it = node_to_unitig_.find(node_id);
    return (it != node_to_unitig_.end()) ? it->second : UINT64_MAX;
}

bool UnitigGraph::detect_cycles() const {
    // DFS-based cycle detection on unitig graph
    size_t n = unitigs_.size();
    enum Color { WHITE, GRAY, BLACK };
    std::vector<Color> color(n, WHITE);

    // Build adjacency
    std::vector<std::vector<uint64_t>> adj(n);
    for (const auto& e : edges_) {
        adj[e.from].push_back(e.to);
    }

    for (size_t start = 0; start < n; ++start) {
        if (color[start] != WHITE) continue;
        std::stack<std::pair<uint64_t, size_t>> stk;
        stk.push({start, 0});
        color[start] = GRAY;
        while (!stk.empty()) {
            auto& [u, idx] = stk.top();
            if (idx < adj[u].size()) {
                uint64_t v = adj[u][idx++];
                if (color[v] == GRAY) return true; // cycle
                if (color[v] == WHITE) {
                    color[v] = GRAY;
                    stk.push({v, 0});
                }
            } else {
                color[u] = BLACK;
                stk.pop();
            }
        }
    }
    return false;
}

bool UnitigGraph::build(const DBG& source) {
    unitigs_.clear();
    edges_.clear();
    hap_edges_.clear();
    node_to_unitig_.clear();

    const auto& nodes = source.nodes();
    size_t num_nodes = nodes.size();
    size_t active_nodes = source.active_node_count();

    // Mark which nodes are "internal" (in-degree 1 AND out-degree 1)
    // and haven't been assigned to a unitig yet.
    std::vector<bool> visited(num_nodes, false);

    // Find maximal non-branching paths
    for (size_t i = 0; i < num_nodes; ++i) {
        const auto& nd = nodes[i];
        if (source.is_node_removed(nd.id)) {
            visited[i] = true;
            continue;
        }
        const auto& in_e  = source.in_edges(nd.id);
        const auto& out_e = source.out_edges(nd.id);

        // Start a unitig from nodes that are NOT internal
        // (i.e., in_degree != 1 or out_degree != 1)
        bool is_internal = (in_e.size() == 1 && out_e.size() == 1);
        if (is_internal || visited[i]) continue;

        // This is a unitig start point. Trace backward and forward along the
        // same linear chain so compaction does not depend on node ID order.
        std::vector<uint64_t> path;
        uint64_t cur = nd.id;

        path.push_back(cur);
        visited[cur] = true;

        // Extend backward first. This handles linear chains whose sink or
        // branch endpoint has a lower node ID than the upstream nodes.
        while (true) {
            const auto& ie = source.in_edges(cur);
            if (ie.size() != 1) break;
            uint64_t prev = source.edges()[ie[0]].from;
            if (visited[prev]) break;
            const auto& prev_out = source.out_edges(prev);
            if (prev_out.size() != 1) break;
            path.insert(path.begin(), prev);
            visited[prev] = true;
            const auto& prev_in = source.in_edges(prev);
            if (prev_in.size() != 1) break;
            cur = prev;
        }

        // Extend forward
        cur = nd.id;
        while (true) {
            const auto& oe = source.out_edges(cur);
            if (oe.size() != 1) break;
            uint64_t next = source.edges()[oe[0]].to;
            if (visited[next]) break;
            const auto& next_in = source.in_edges(next);
            if (next_in.size() != 1) break; // next is a branch point, don't consume
            // Also stop if next has out-degree != 1 (it's an endpoint)
            const auto& next_out = source.out_edges(next);
            path.push_back(next);
            visited[next] = true;
            if (next_out.size() != 1) break;
            cur = next;
        }

        // Create unitig
        Unitig u;
        u.id = unitigs_.size();
        u.node_ids = path;

        // Build sequence: first node's kmer + last char of each subsequent node's kmer
        u.sequence = source.node(path[0]).kmer;
        for (size_t j = 1; j < path.size(); ++j) {
            const auto& km = source.node(path[j]).kmer;
            u.sequence += km.back();
        }

        // Mean depth
        double total = 0;
        for (uint64_t nid : path) total += source.node(nid).depth;
        u.mean_depth = total / path.size();

        // Map nodes to unitig
        for (uint64_t nid : path) {
            node_to_unitig_[nid] = u.id;
        }

        unitigs_.push_back(std::move(u));
    }

    // Also handle isolated internal-only chains that were missed
    // (all-internal cycles would be caught here, but we detect cycles later)
    for (size_t i = 0; i < num_nodes; ++i) {
        if (visited[i]) continue;
        const auto& nd = nodes[i];
        if (source.is_node_removed(nd.id)) continue;
        std::vector<uint64_t> path;
        uint64_t cur = nd.id;

        while (!visited[cur]) {
            path.push_back(cur);
            visited[cur] = true;
            const auto& oe = source.out_edges(cur);
            if (oe.size() != 1) break;
            cur = source.edges()[oe[0]].to;
        }

        if (path.empty()) continue;

        Unitig u;
        u.id = unitigs_.size();
        u.node_ids = path;
        u.sequence = source.node(path[0]).kmer;
        for (size_t j = 1; j < path.size(); ++j)
            u.sequence += source.node(path[j]).kmer.back();
        double total = 0;
        for (uint64_t nid : path) total += source.node(nid).depth;
        u.mean_depth = total / path.size();
        for (uint64_t nid : path) node_to_unitig_[nid] = u.id;
        unitigs_.push_back(std::move(u));
    }

    std::vector<uint64_t> predecessor(unitigs_.size(), UINT64_MAX);
    std::vector<uint64_t> successor(unitigs_.size(), UINT64_MAX);
    for (size_t unitig_id = 0; unitig_id < unitigs_.size(); ++unitig_id) {
        const auto& unitig = unitigs_[unitig_id];

        const uint64_t first_node = unitig.node_ids.front();
        const auto& in_edges = source.in_edges(first_node);
        if (in_edges.size() == 1) {
            const uint64_t prev_node = source.edges()[in_edges[0]].from;
            const uint64_t prev_unitig = node_to_unitig(prev_node);
            if (prev_unitig != UINT64_MAX && prev_unitig != unitig_id &&
                source.out_edges(prev_node).size() == 1) {
                predecessor[unitig_id] = prev_unitig;
            }
        }

        const uint64_t last_node = unitig.node_ids.back();
        const auto& out_edges = source.out_edges(last_node);
        if (out_edges.size() == 1) {
            const uint64_t next_node = source.edges()[out_edges[0]].to;
            const uint64_t next_unitig = node_to_unitig(next_node);
            if (next_unitig != UINT64_MAX && next_unitig != unitig_id &&
                source.in_edges(next_node).size() == 1) {
                successor[unitig_id] = next_unitig;
            }
        }
    }

    auto build_unitig_from_path = [&source](uint64_t unitig_id,
                                            const std::vector<uint64_t>& path) {
        Unitig unitig;
        unitig.id = unitig_id;
        unitig.node_ids = path;
        unitig.sequence = source.node(path[0]).kmer;
        for (size_t index = 1; index < path.size(); ++index) {
            unitig.sequence += source.node(path[index]).kmer.back();
        }

        double total_depth = 0.0;
        for (uint64_t node_id : path) {
            total_depth += source.node(node_id).depth;
        }
        unitig.mean_depth = total_depth / path.size();
        return unitig;
    };

    std::vector<bool> merged_visited(unitigs_.size(), false);
    std::vector<Unitig> merged_unitigs;
    merged_unitigs.reserve(unitigs_.size());

    auto append_chain = [&](uint64_t start_unitig) {
        std::vector<uint64_t> merged_path;
        uint64_t current = start_unitig;

        while (!merged_visited[current]) {
            merged_visited[current] = true;
            const auto& current_nodes = unitigs_[current].node_ids;
            merged_path.insert(merged_path.end(), current_nodes.begin(), current_nodes.end());

            const uint64_t next = successor[current];
            if (next == UINT64_MAX || predecessor[next] != current) {
                break;
            }
            current = next;
        }

        merged_unitigs.push_back(build_unitig_from_path(merged_unitigs.size(), merged_path));
    };

    for (size_t unitig_id = 0; unitig_id < unitigs_.size(); ++unitig_id) {
        if (merged_visited[unitig_id] || predecessor[unitig_id] != UINT64_MAX) {
            continue;
        }
        append_chain(unitig_id);
    }

    for (size_t unitig_id = 0; unitig_id < unitigs_.size(); ++unitig_id) {
        if (merged_visited[unitig_id]) {
            continue;
        }
        merged_visited[unitig_id] = true;
        merged_unitigs.push_back(
            build_unitig_from_path(merged_unitigs.size(), unitigs_[unitig_id].node_ids));
    }

    if (merged_unitigs.size() != unitigs_.size()) {
        unitigs_ = std::move(merged_unitigs);
        node_to_unitig_.clear();
        for (const auto& unitig : unitigs_) {
            for (uint64_t node_id : unitig.node_ids) {
                node_to_unitig_[node_id] = unitig.id;
            }
        }
    }

    spdlog::info("Unitig graph: {} unitigs from {} active nodes ({} total nodes)",
                 unitigs_.size(), active_nodes, num_nodes);

    // Build unitig-level edges from the original edges
    // An edge between unitigs exists when the last node of one unitig
    // connects to the first node of another (or internal connections across unitig boundaries).
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, uint32_t>> ug_edge_map;
    for (const auto& e : source.edges()) {
        uint64_t u_from = node_to_unitig(e.from);
        uint64_t u_to   = node_to_unitig(e.to);
        if (u_from == UINT64_MAX || u_to == UINT64_MAX) continue;
        if (u_from == u_to) continue; // internal to same unitig
        ug_edge_map[u_from][u_to] += e.weight;
    }
    for (const auto& [from, targets] : ug_edge_map) {
        for (const auto& [to, weight] : targets) {
            edges_.push_back({from, to, weight});
        }
    }

    // Build unitig-level haplotype edges
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, uint32_t>> ug_hap_map;
    for (const auto& he : source.haplotype_edges()) {
        uint64_t u_from = node_to_unitig(he.from_node);
        uint64_t u_to   = node_to_unitig(he.to_node);
        if (u_from == UINT64_MAX || u_to == UINT64_MAX) continue;
        if (u_from == u_to) continue;
        ug_hap_map[u_from][u_to] += he.weight;
    }
    for (const auto& [from, targets] : ug_hap_map) {
        for (const auto& [to, weight] : targets) {
            hap_edges_.push_back({from, to, weight});
        }
    }

    // Prune low-weight haplotype edges: weight < 5% of min(cov_U1, cov_U2)
    hap_edges_.erase(
        std::remove_if(hap_edges_.begin(), hap_edges_.end(),
            [this](const HaplotypeEdge& he) {
                double min_cov = std::min(unitigs_[he.from_node].mean_depth,
                                          unitigs_[he.to_node].mean_depth);
                return he.weight < 0.05 * min_cov;
            }),
        hap_edges_.end());

    std::vector<size_t> graph_incident_counts(unitigs_.size(), 0);
    for (const auto& edge : edges_) {
        graph_incident_counts[edge.from]++;
        graph_incident_counts[edge.to]++;
    }

    std::vector<uint64_t> unitig_id_remap(unitigs_.size(), UINT64_MAX);
    std::vector<Unitig> filtered_unitigs;
    filtered_unitigs.reserve(unitigs_.size());
    for (size_t old_id = 0; old_id < unitigs_.size(); ++old_id) {
        const auto& unitig = unitigs_[old_id];
        const bool is_isolated_singleton =
            unitig.node_ids.size() == 1 && graph_incident_counts[old_id] == 0;
        if (is_isolated_singleton) {
            continue;
        }

        Unitig kept = unitig;
        kept.id = filtered_unitigs.size();
        unitig_id_remap[old_id] = kept.id;
        filtered_unitigs.push_back(std::move(kept));
    }

    if (filtered_unitigs.size() != unitigs_.size()) {
        std::vector<Edge> remapped_edges;
        remapped_edges.reserve(edges_.size());
        for (const auto& edge : edges_) {
            if (unitig_id_remap[edge.from] == UINT64_MAX || unitig_id_remap[edge.to] == UINT64_MAX) {
                continue;
            }
            remapped_edges.push_back({unitig_id_remap[edge.from], unitig_id_remap[edge.to], edge.weight});
        }

        std::vector<HaplotypeEdge> remapped_hap_edges;
        remapped_hap_edges.reserve(hap_edges_.size());
        for (const auto& hap_edge : hap_edges_) {
            if (unitig_id_remap[hap_edge.from_node] == UINT64_MAX ||
                unitig_id_remap[hap_edge.to_node] == UINT64_MAX) {
                continue;
            }
            remapped_hap_edges.push_back(
                {unitig_id_remap[hap_edge.from_node], unitig_id_remap[hap_edge.to_node], hap_edge.weight});
        }

        unitigs_ = std::move(filtered_unitigs);
        edges_ = std::move(remapped_edges);
        hap_edges_ = std::move(remapped_hap_edges);
        node_to_unitig_.clear();
        for (const auto& unitig : unitigs_) {
            for (uint64_t node_id : unitig.node_ids) {
                node_to_unitig_[node_id] = unitig.id;
            }
        }
    }

    spdlog::info("Unitig graph: {} edges, {} haplotype edges",
                 edges_.size(), hap_edges_.size());

    // Cycle detection
    has_cycles_ = detect_cycles();
    if (has_cycles_) {
        spdlog::warn("Cycle detected in unitig graph — aborting assembly");
        return false;
    }

    return true;
}

} // namespace sharda
