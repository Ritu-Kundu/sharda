#include "assembly/graph_cleaner.h"
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace sharda {

namespace {

/// Follow a dead-end path from `start_node` in direction `forward`.
/// Returns the list of node IDs on the tip (excluding branch point).
std::vector<uint64_t> trace_tip(const DBG& graph, uint64_t start_node, bool forward) {
    std::vector<uint64_t> tip;
    uint64_t cur = start_node;
    while (true) {
        tip.push_back(cur);
        const auto& adj = forward ? graph.out_edges(cur) : graph.in_edges(cur);
        if (adj.size() != 1) break;
        uint64_t next_edge = adj[0];
        const auto& e = graph.edges()[next_edge];
        uint64_t next = forward ? e.to : e.from;
        // Check the other direction: if next has branching, stop
        const auto& rev = forward ? graph.in_edges(next) : graph.out_edges(next);
        if (rev.size() > 1) break; // next is a branch point, don't include
        cur = next;
    }
    return tip;
}

bool tip_has_haplotype_support(const DBG& graph, const std::vector<uint64_t>& tip) {
    if (tip.empty()) return false;

    std::unordered_set<uint64_t> tip_nodes(tip.begin(), tip.end());
    for (const auto& he : graph.haplotype_edges()) {
        if (tip_nodes.count(he.from_node) || tip_nodes.count(he.to_node)) {
            return true;
        }
    }
    return false;
}

/// Remove tips: dead-end paths shorter than threshold.
/// Returns number of nodes removed.
int remove_tips(DBG& graph, int max_tip_len) {
    int removed = 0;
    for (const auto& node : graph.nodes()) {
        uint64_t nid = node.id;
        if (node.is_backbone) continue;

        // Check for dead-end at start (in-degree == 0, out-degree > 0)
        if (graph.in_edges(nid).empty() && !graph.out_edges(nid).empty()) {
            auto tip = trace_tip(graph, nid, true);
            if (static_cast<int>(tip.size()) < max_tip_len &&
                !tip_has_haplotype_support(graph, tip)) {
                for (uint64_t id : tip) graph.remove_node(id);
                removed += tip.size();
            }
        }
        // Check for dead-end at end (out-degree == 0, in-degree > 0)
        if (graph.out_edges(nid).empty() && !graph.in_edges(nid).empty()) {
            auto tip = trace_tip(graph, nid, false);
            if (static_cast<int>(tip.size()) < max_tip_len &&
                !tip_has_haplotype_support(graph, tip)) {
                for (uint64_t id : tip) graph.remove_node(id);
                removed += tip.size();
            }
        }
    }
    return removed;
}

/// Compute local average edge weight around a reference position (±window).
double local_avg_weight(const DBG& graph, int32_t ref_pos, int window = 500) {
    double sum = 0;
    int count = 0;
    for (const auto& e : graph.edges()) {
        const auto& from_node = graph.node(e.from);
        int32_t pos = from_node.ref_pos;
        if (pos >= 0 && std::abs(pos - ref_pos) <= window) {
            sum += e.weight;
            count++;
        }
    }
    return (count > 0) ? sum / count : 1.0;
}

/// Prune low-weight edges (< 5% of local average).
/// Returns number of edges removed.
int prune_low_weight_edges(DBG& graph) {
    int removed = 0;
    for (size_t i = 0; i < graph.edges().size(); ++i) {
        const auto& e = graph.edges()[i];
        const auto& from_node = graph.node(e.from);
        const auto& to_node = graph.node(e.to);

        // Preserve read-supported branch edges. In low-coverage examples,
        // variant paths often carry single-read support against a higher-depth
        // backbone, so pruning them by local backbone coverage erases the only
        // alternative path before compaction.
        if (!from_node.is_backbone || !to_node.is_backbone) {
            continue;
        }

        int32_t pos = from_node.ref_pos;
        if (pos < 0) {
            // For non-backbone, use the to-node's pos
            pos = to_node.ref_pos;
        }
        double avg = local_avg_weight(graph, pos);
        if (e.weight < 0.05 * avg) {
            graph.remove_edge(i);
            removed++;
        }
    }
    return removed;
}

/// Simple bubble popping: find diverge-reconverge patterns and remove weak branch.
/// Returns number of nodes removed.
int pop_bubbles(DBG& graph) {
    int removed = 0;
    int k = graph.k();

    for (const auto& node : graph.nodes()) {
        const auto& out = graph.out_edges(node.id);
        if (out.size() != 2) continue;

        // Two branches diverge from this node
        const auto& e0 = graph.edges()[out[0]];
        const auto& e1 = graph.edges()[out[1]];
        uint64_t b0 = e0.to, b1 = e1.to;

        // Trace each branch forward to see if they reconverge within 2*k
        auto trace_forward = [&](uint64_t start) -> std::vector<uint64_t> {
            std::vector<uint64_t> path;
            uint64_t cur = start;
            for (int step = 0; step < 2 * k; ++step) {
                path.push_back(cur);
                const auto& oe = graph.out_edges(cur);
                if (oe.size() != 1) break;
                cur = graph.edges()[oe[0]].to;
            }
            return path;
        };

        auto path0 = trace_forward(b0);
        auto path1 = trace_forward(b1);

        if (path0.empty() || path1.empty()) continue;

        uint64_t end0 = path0.back();
        uint64_t end1 = path1.back();

        // Check reconvergence: do they end at the same node or share a successor?
        if (end0 == end1 || (!graph.out_edges(end0).empty() && !graph.out_edges(end1).empty() &&
            graph.edges()[graph.out_edges(end0)[0]].to == graph.edges()[graph.out_edges(end1)[0]].to)) {
            // Bubble detected. Compute weights of branches.
            uint32_t w0 = e0.weight;
            uint32_t w1 = e1.weight;

            // Remove the weaker branch if weight < 20% of the stronger
            auto remove_path = [&](const std::vector<uint64_t>& path) {
                for (uint64_t nid : path) {
                    if (!graph.node(nid).is_backbone)
                        graph.remove_node(nid);
                }
                removed += path.size();
            };

            if (w0 < w1 && w0 < 0.2 * w1 && path0.size() <= path1.size() + 2) {
                remove_path(path0);
            } else if (w1 < w0 && w1 < 0.2 * w0 && path1.size() <= path0.size() + 2) {
                remove_path(path1);
            }
        }
    }
    return removed;
}

} // anonymous namespace

void clean_graph(DBG& graph, int mean_read_length) {
    spdlog::info("Graph cleaning started: {} nodes, {} edges",
                 graph.node_count(), graph.edge_count());

    for (int iteration = 0; iteration < 10; ++iteration) {
        size_t edges_before_tips = graph.edge_count();
        int tips     = remove_tips(graph, mean_read_length);
        graph.rebuild_adjacency();
        size_t edges_after_tips = graph.edge_count();
        size_t tip_edges_removed = edges_before_tips - edges_after_tips;

        size_t edges_before_low_wt = graph.edge_count();
        int low_wt   = prune_low_weight_edges(graph);
        graph.rebuild_adjacency();
        size_t edges_after_low_wt = graph.edge_count();
        size_t low_wt_edges_removed = edges_before_low_wt - edges_after_low_wt;

        size_t edges_before_bubbles = graph.edge_count();
        int bubbles  = pop_bubbles(graph);
        graph.rebuild_adjacency();
        size_t edges_after_bubbles = graph.edge_count();
        size_t bubble_edges_removed = edges_before_bubbles - edges_after_bubbles;

        spdlog::info(
            "Cleaning iteration {}: tip_nodes_removed={}, tip_edges_removed={}, "
            "low_weight_edges_flagged={}, low_weight_edges_removed={}, "
            "bubble_nodes_removed={}, bubble_edges_removed={}, remaining_edges={}",
            iteration, tips, tip_edges_removed,
            low_wt, low_wt_edges_removed,
            bubbles, bubble_edges_removed,
            graph.edge_count());

        if (tips == 0 && low_wt == 0 && bubbles == 0) break;
    }

    spdlog::info("Graph cleaning done: {} nodes, {} edges",
                 graph.node_count(), graph.edge_count());
}

} // namespace sharda
