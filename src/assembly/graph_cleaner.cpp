#include "assembly/graph_cleaner.h"
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
#include <sstream>

namespace sharda {

namespace {

constexpr double kRelativeSupportThreshold = 0.05;
constexpr double kRegionMeanDepthThreshold = 0.25;
constexpr int kLocalCoverageWindow = 500;

std::string format_tip_nodes(const std::vector<uint64_t>& tip) {
    std::ostringstream out;
    out << '[';
    for (size_t index = 0; index < tip.size(); ++index) {
        if (index > 0) {
            out << ", ";
        }
        out << tip[index];
    }
    out << ']';
    return out.str();
}

std::vector<int32_t> collect_node_positions(const Node& node) {
    if (!node.ref_positions.empty()) {
        return {node.ref_positions.begin(), node.ref_positions.end()};
    }

    if (node.ref_pos >= 0) {
        return {node.ref_pos};
    }

    return {};
}

std::vector<int32_t> collect_edge_positions(const DBG& graph, const Edge& edge) {
    std::vector<int32_t> positions = collect_node_positions(graph.node(edge.from));
    const auto to_positions = collect_node_positions(graph.node(edge.to));

    for (int32_t pos : to_positions) {
        if (std::find(positions.begin(), positions.end(), pos) == positions.end()) {
            positions.push_back(pos);
        }
    }

    return positions;
}

bool positions_overlap_window(const std::vector<int32_t>& lhs,
                              const std::vector<int32_t>& rhs,
                              int window) {
    for (int32_t left : lhs) {
        for (int32_t right : rhs) {
            if (std::abs(left - right) <= window) {
                return true;
            }
        }
    }
    return false;
}

bool edge_within_window(const DBG& graph,
                        const Edge& edge,
                        const std::vector<int32_t>& anchor_positions,
                        int window) {
    if (anchor_positions.empty()) {
        return false;
    }

    const auto edge_positions = collect_edge_positions(graph, edge);
    return positions_overlap_window(edge_positions, anchor_positions, window);
}

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

double local_avg_weight(const DBG& graph,
                        const std::vector<int32_t>& anchor_positions,
                        int window = kLocalCoverageWindow) {
    if (anchor_positions.empty()) {
        return 1.0;
    }

    double sum = 0;
    int count = 0;
    for (const auto& edge : graph.edges()) {
        if (edge_within_window(graph, edge, anchor_positions, window)) {
            sum += edge.weight;
            count++;
        }
    }

    return (count > 0) ? sum / count : 1.0;
}

double mean_backbone_depth(const DBG& graph) {
    double sum = 0.0;
    int count = 0;
    for (const auto& node : graph.nodes()) {
        if (!node.is_backbone || graph.is_node_removed(node.id)) {
            continue;
        }
        sum += node.depth;
        count++;
    }

    return (count > 0) ? sum / count : 0.0;
}

double support_threshold(const DBG& graph, double local_avg, bool use_regional_floor) {
    const double relative_threshold = kRelativeSupportThreshold * local_avg;
    if (!use_regional_floor) {
        return relative_threshold;
    }

    const double regional_floor = mean_backbone_depth(graph) * kRegionMeanDepthThreshold;
    return std::max(relative_threshold, regional_floor);
}

uint32_t edge_weight_between(const DBG& graph, uint64_t from, uint64_t to) {
    for (uint64_t edge_index : graph.out_edges(from)) {
        const auto& edge = graph.edges()[edge_index];
        if (edge.to == to) {
            return edge.weight;
        }
    }
    return 0;
}

std::vector<uint32_t> collect_tip_edge_weights(const DBG& graph,
                                               const std::vector<uint64_t>& tip,
                                               bool forward) {
    std::vector<uint32_t> weights;
    if (tip.empty()) {
        return weights;
    }

    for (size_t i = 0; i + 1 < tip.size(); ++i) {
        uint64_t from = forward ? tip[i] : tip[i + 1];
        uint64_t to = forward ? tip[i + 1] : tip[i];
        uint32_t weight = edge_weight_between(graph, from, to);
        if (weight > 0) {
            weights.push_back(weight);
        }
    }

    uint64_t boundary_node = tip.back();
    const auto& boundary_edges = forward ? graph.out_edges(boundary_node)
                                         : graph.in_edges(boundary_node);
    if (boundary_edges.size() == 1) {
        weights.push_back(graph.edges()[boundary_edges[0]].weight);
    }

    return weights;
}

std::vector<int32_t> collect_tip_positions(const DBG& graph,
                                           const std::vector<uint64_t>& tip) {
    std::vector<int32_t> positions;
    for (uint64_t node_id : tip) {
        const auto node_positions = collect_node_positions(graph.node(node_id));
        for (int32_t pos : node_positions) {
            if (std::find(positions.begin(), positions.end(), pos) == positions.end()) {
                positions.push_back(pos);
            }
        }
    }
    return positions;
}

double mean_weight(const std::vector<uint32_t>& weights) {
    if (weights.empty()) {
        return 0.0;
    }

    double sum = 0.0;
    for (uint32_t weight : weights) {
        sum += weight;
    }
    return sum / weights.size();
}

struct TipMetrics {
    double tip_support = 0.0;
    double local_avg = 1.0;
};

TipMetrics compute_tip_metrics(const DBG& graph,
                               const std::vector<uint64_t>& tip,
                               bool forward) {
    TipMetrics metrics;
    metrics.local_avg = local_avg_weight(graph, collect_tip_positions(graph, tip));
    metrics.tip_support = mean_weight(collect_tip_edge_weights(graph, tip, forward));
    return metrics;
}

void log_tip_rejection(uint64_t start_node,
                       bool forward,
                       const std::vector<uint64_t>& tip,
                       const char* reason,
                       const TipMetrics& metrics,
                       double threshold,
                       int max_tip_len) {
    spdlog::debug(
        "Tip prune rejected: start_node={} direction={} chain={} length={} max_tip_len={} "
        "tip_support={} local_avg={} threshold={} reason={}",
        start_node,
        forward ? "forward" : "reverse",
        format_tip_nodes(tip),
        tip.size(),
        max_tip_len,
        metrics.tip_support,
        metrics.local_avg,
        threshold,
        reason);
}

void log_non_tip_skip(const DBG& graph, const Node& node) {
    const auto in_degree = graph.in_edges(node.id).size();
    const auto out_degree = graph.out_edges(node.id).size();
    if (node.depth > 1 || (in_degree > 1 && out_degree > 1)) {
        return;
    }

    spdlog::debug(
        "Tip prune skipped node {}: not a dead-end (in_degree={}, out_degree={}, depth={})",
        node.id,
        in_degree,
        out_degree,
        node.depth);
}

/// Remove tips: dead-end paths shorter than threshold and under-supported
/// relative to local coverage. Returns number of nodes removed.
int remove_tips(DBG& graph, int max_tip_len) {
    int removed = 0;
    for (const auto& node : graph.nodes()) {
        uint64_t nid = node.id;
        if (node.is_backbone || graph.is_node_removed(nid)) continue;

        const bool is_forward_tip = graph.in_edges(nid).empty() && !graph.out_edges(nid).empty();
        const bool is_reverse_tip = graph.out_edges(nid).empty() && !graph.in_edges(nid).empty();

        // Check for dead-end at start (in-degree == 0, out-degree > 0)
        if (is_forward_tip) {
            auto tip = trace_tip(graph, nid, true);
            TipMetrics metrics = compute_tip_metrics(graph, tip, true);
            const double threshold = support_threshold(graph, metrics.local_avg, true);
            if (static_cast<int>(tip.size()) < max_tip_len &&
                metrics.tip_support < threshold) {
                for (uint64_t id : tip) {
                    if (!graph.is_node_removed(id)) {
                        graph.remove_node(id);
                        removed++;
                    }
                }
                spdlog::debug(
                    "Tip prune removed: start_node={} direction=forward chain={} length={} "
                    "tip_support={} local_avg={} threshold={}",
                    nid,
                    format_tip_nodes(tip),
                    tip.size(),
                    metrics.tip_support,
                    metrics.local_avg,
                    threshold);
            } else if (static_cast<int>(tip.size()) >= max_tip_len) {
                log_tip_rejection(nid, true, tip, "length", metrics, threshold, max_tip_len);
            } else {
                log_tip_rejection(nid, true, tip, "support", metrics, threshold, max_tip_len);
            }
        }
        // Check for dead-end at end (out-degree == 0, in-degree > 0)
        if (is_reverse_tip) {
            auto tip = trace_tip(graph, nid, false);
            TipMetrics metrics = compute_tip_metrics(graph, tip, false);
            const double threshold = support_threshold(graph, metrics.local_avg, true);
            if (static_cast<int>(tip.size()) < max_tip_len &&
                metrics.tip_support < threshold) {
                for (uint64_t id : tip) {
                    if (!graph.is_node_removed(id)) {
                        graph.remove_node(id);
                        removed++;
                    }
                }
                spdlog::debug(
                    "Tip prune removed: start_node={} direction=reverse chain={} length={} "
                    "tip_support={} local_avg={} threshold={}",
                    nid,
                    format_tip_nodes(tip),
                    tip.size(),
                    metrics.tip_support,
                    metrics.local_avg,
                    threshold);
            } else if (static_cast<int>(tip.size()) >= max_tip_len) {
                log_tip_rejection(nid, false, tip, "length", metrics, threshold, max_tip_len);
            } else {
                log_tip_rejection(nid, false, tip, "support", metrics, threshold, max_tip_len);
            }
        }

        if (!is_forward_tip && !is_reverse_tip) {
            log_non_tip_skip(graph, node);
        }
    }
    return removed;
}

/// Prune low-weight edges (< 5% of local average).
/// Returns number of edges removed.
int prune_low_weight_edges(DBG& graph) {
    int removed = 0;
    for (size_t i = 0; i < graph.edges().size(); ++i) {
        const auto& e = graph.edges()[i];
        const auto anchor_positions = collect_edge_positions(graph, e);
        double avg = local_avg_weight(graph, anchor_positions);
        const bool uses_non_backbone = !graph.node(e.from).is_backbone || !graph.node(e.to).is_backbone;
        if (e.weight < support_threshold(graph, avg, uses_non_backbone)) {
            graph.remove_edge(i);
            removed++;
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

        spdlog::info(
            "Cleaning iteration {}: tip_nodes_removed={}, tip_edges_removed={}, "
            "low_weight_edges_flagged={}, low_weight_edges_removed={}, "
            "bubble_popping_skipped=true, remaining_edges={}",
            iteration, tips, tip_edges_removed,
            low_wt, low_wt_edges_removed,
            graph.edge_count());

        if (tips == 0 && low_wt == 0) break;
    }

    spdlog::info("Graph cleaning done: {} nodes, {} edges",
                 graph.node_count(), graph.edge_count());
}

} // namespace sharda
