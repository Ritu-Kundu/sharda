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
constexpr int kMinInternalBranchNodes = 2;

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

std::vector<uint32_t> collect_path_edge_weights(const DBG& graph,
                                                const std::vector<uint64_t>& path,
                                                uint64_t source_anchor,
                                                uint64_t sink_anchor) {
    std::vector<uint32_t> weights;
    if (path.empty()) {
        return weights;
    }

    if (source_anchor != UINT64_MAX) {
        uint32_t weight = edge_weight_between(graph, source_anchor, path.front());
        if (weight > 0) {
            weights.push_back(weight);
        }
    }

    for (size_t i = 0; i + 1 < path.size(); ++i) {
        uint32_t weight = edge_weight_between(graph, path[i], path[i + 1]);
        if (weight > 0) {
            weights.push_back(weight);
        }
    }

    if (sink_anchor != UINT64_MAX) {
        uint32_t weight = edge_weight_between(graph, path.back(), sink_anchor);
        if (weight > 0) {
            weights.push_back(weight);
        }
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

bool is_internal_non_backbone_node(const DBG& graph, uint64_t node_id) {
    const auto& node = graph.node(node_id);
    return !node.is_backbone &&
           !graph.is_node_removed(node_id) &&
           graph.in_edges(node_id).size() == 1 &&
           graph.out_edges(node_id).size() == 1;
}

bool is_active_non_backbone_node(const DBG& graph, uint64_t node_id) {
    const auto& node = graph.node(node_id);
    return !node.is_backbone && !graph.is_node_removed(node_id);
}

struct InternalBranch {
    std::vector<uint64_t> nodes;
    uint64_t source_anchor = UINT64_MAX;
    uint64_t sink_anchor = UINT64_MAX;
    double branch_support = 0.0;
    double local_avg = 1.0;
};

InternalBranch trace_internal_branch(const DBG& graph, uint64_t start_node) {
    InternalBranch branch;
    if (!is_internal_non_backbone_node(graph, start_node)) {
        return branch;
    }

    branch.nodes.push_back(start_node);

    while (true) {
        const auto& in_edges = graph.in_edges(branch.nodes.front());
        if (in_edges.size() != 1) {
            break;
        }
        const uint64_t prev = graph.edges()[in_edges[0]].from;
        if (!is_internal_non_backbone_node(graph, prev)) {
            break;
        }
        branch.nodes.insert(branch.nodes.begin(), prev);
    }

    while (true) {
        const auto& out_edges = graph.out_edges(branch.nodes.back());
        if (out_edges.size() != 1) {
            break;
        }
        const uint64_t next = graph.edges()[out_edges[0]].to;
        if (!is_internal_non_backbone_node(graph, next)) {
            break;
        }
        branch.nodes.push_back(next);
    }

    const auto& source_in_edges = graph.in_edges(branch.nodes.front());
    const auto& sink_out_edges = graph.out_edges(branch.nodes.back());
    if (source_in_edges.size() != 1 || sink_out_edges.size() != 1) {
        branch.nodes.clear();
        return branch;
    }

    branch.source_anchor = graph.edges()[source_in_edges[0]].from;
    branch.sink_anchor = graph.edges()[sink_out_edges[0]].to;
    if (branch.source_anchor == UINT64_MAX || branch.sink_anchor == UINT64_MAX) {
        branch.nodes.clear();
        return branch;
    }

    const bool has_divergence = graph.out_edges(branch.source_anchor).size() > 1;
    const bool has_convergence = graph.in_edges(branch.sink_anchor).size() > 1;
    if (!has_divergence || !has_convergence) {
        branch.nodes.clear();
        return branch;
    }

    branch.local_avg = local_avg_weight(graph, collect_tip_positions(graph, branch.nodes));
    branch.branch_support = mean_weight(
        collect_path_edge_weights(graph, branch.nodes, branch.source_anchor, branch.sink_anchor));
    return branch;
}

InternalBranch trace_internal_component(const DBG& graph, uint64_t start_node) {
    InternalBranch branch;
    if (!is_active_non_backbone_node(graph, start_node)) {
        return branch;
    }

    std::vector<bool> in_component(graph.node_count(), false);
    std::vector<uint64_t> stack = {start_node};
    std::vector<uint64_t> source_candidates;
    std::vector<uint64_t> sink_candidates;

    while (!stack.empty()) {
        const uint64_t node_id = stack.back();
        stack.pop_back();
        if (node_id >= in_component.size() || in_component[node_id] || !is_active_non_backbone_node(graph, node_id)) {
            continue;
        }

        in_component[node_id] = true;
        branch.nodes.push_back(node_id);

        for (uint64_t edge_index : graph.in_edges(node_id)) {
            const auto& edge = graph.edges()[edge_index];
            if (edge.weight > 1) {
                branch.nodes.clear();
                return branch;
            }

            if (is_active_non_backbone_node(graph, edge.from)) {
                stack.push_back(edge.from);
            } else {
                source_candidates.push_back(edge.from);
            }
        }

        for (uint64_t edge_index : graph.out_edges(node_id)) {
            const auto& edge = graph.edges()[edge_index];
            if (edge.weight > 1) {
                branch.nodes.clear();
                return branch;
            }

            if (is_active_non_backbone_node(graph, edge.to)) {
                stack.push_back(edge.to);
            } else {
                sink_candidates.push_back(edge.to);
            }
        }
    }

    if (branch.nodes.empty()) {
        return branch;
    }

    std::sort(branch.nodes.begin(), branch.nodes.end());
    branch.nodes.erase(std::unique(branch.nodes.begin(), branch.nodes.end()), branch.nodes.end());
    std::sort(source_candidates.begin(), source_candidates.end());
    source_candidates.erase(std::unique(source_candidates.begin(), source_candidates.end()), source_candidates.end());
    std::sort(sink_candidates.begin(), sink_candidates.end());
    sink_candidates.erase(std::unique(sink_candidates.begin(), sink_candidates.end()), sink_candidates.end());

    const bool has_divergence_anchor = std::any_of(
        source_candidates.begin(), source_candidates.end(), [&graph](uint64_t node_id) {
            return !graph.is_node_removed(node_id) && graph.out_edges(node_id).size() > 1;
        });
    const bool has_convergence_anchor = std::any_of(
        sink_candidates.begin(), sink_candidates.end(), [&graph](uint64_t node_id) {
            return !graph.is_node_removed(node_id) && graph.in_edges(node_id).size() > 1;
        });
    if (!has_divergence_anchor && !has_convergence_anchor) {
        branch.nodes.clear();
        return branch;
    }

    branch.source_anchor = source_candidates.empty() ? UINT64_MAX : source_candidates.front();
    branch.sink_anchor = sink_candidates.empty() ? UINT64_MAX : sink_candidates.front();

    std::vector<int32_t> positions;
    for (uint64_t node_id : branch.nodes) {
        const auto node_positions = collect_node_positions(graph.node(node_id));
        for (int32_t pos : node_positions) {
            if (std::find(positions.begin(), positions.end(), pos) == positions.end()) {
                positions.push_back(pos);
            }
        }
    }
    branch.local_avg = local_avg_weight(graph, positions);

    std::vector<uint32_t> weights;
    for (uint64_t node_id : branch.nodes) {
        for (uint64_t edge_index : graph.in_edges(node_id)) {
            const auto& edge = graph.edges()[edge_index];
            if ((edge.from < in_component.size() && in_component[edge.from]) || edge.weight == 0) {
                continue;
            }
            weights.push_back(edge.weight);
        }
        for (uint64_t edge_index : graph.out_edges(node_id)) {
            const auto& edge = graph.edges()[edge_index];
            if ((edge.to < in_component.size() && in_component[edge.to]) || edge.weight == 0) {
                continue;
            }
            weights.push_back(edge.weight);
        }
    }
    branch.branch_support = mean_weight(weights);
    return branch;
}

int remove_weak_internal_branches(DBG& graph, int max_branch_len) {
    int removed = 0;
    std::vector<bool> visited(graph.node_count(), false);

    for (const auto& node : graph.nodes()) {
        const uint64_t node_id = node.id;
        if (node_id >= visited.size() || visited[node_id] || !is_active_non_backbone_node(graph, node_id)) {
            continue;
        }

        InternalBranch branch = trace_internal_branch(graph, node_id);
        if (branch.nodes.size() < static_cast<size_t>(kMinInternalBranchNodes)) {
            branch = trace_internal_component(graph, node_id);
        }
        for (uint64_t branch_node : branch.nodes) {
            if (branch_node < visited.size()) {
                visited[branch_node] = true;
            }
        }

        if (branch.nodes.size() < static_cast<size_t>(kMinInternalBranchNodes)) {
            continue;
        }

        const double contextual_threshold = support_threshold(graph, branch.local_avg, true);
        const double threshold = 1.0;
        if (static_cast<int>(branch.nodes.size()) >= max_branch_len ||
            branch.branch_support > threshold ||
            contextual_threshold > threshold) {
            spdlog::debug(
                "Internal branch kept: source_anchor={} sink_anchor={} chain={} length={} "
                "branch_support={} local_avg={} threshold={} contextual_threshold={}",
                branch.source_anchor,
                branch.sink_anchor,
                format_tip_nodes(branch.nodes),
                branch.nodes.size(),
                branch.branch_support,
                branch.local_avg,
                threshold,
                contextual_threshold);
            continue;
        }

        for (uint64_t branch_node : branch.nodes) {
            if (!graph.is_node_removed(branch_node)) {
                graph.remove_node(branch_node);
                removed++;
            }
        }

        spdlog::debug(
            "Internal branch pruned: source_anchor={} sink_anchor={} chain={} length={} "
            "branch_support={} local_avg={} threshold={} contextual_threshold={}",
            branch.source_anchor,
            branch.sink_anchor,
            format_tip_nodes(branch.nodes),
            branch.nodes.size(),
            branch.branch_support,
            branch.local_avg,
            threshold,
            contextual_threshold);
    }

    return removed;
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
int prune_low_weight_edges(DBG& graph, GraphCleaningOptions options) {
    int removed = 0;
    for (size_t i = 0; i < graph.edges().size(); ++i) {
        const auto& e = graph.edges()[i];
        const bool is_backbone_edge = graph.node(e.from).is_backbone && graph.node(e.to).is_backbone;
        if (options.preserve_backbone_edges && is_backbone_edge) {
            continue;
        }
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

void clean_graph(DBG& graph,
                 int mean_read_length,
                 GraphCleaningOptions options) {
    spdlog::info("Graph cleaning started: {} nodes, {} edges",
                 graph.node_count(), graph.edge_count());

    for (int iteration = 0; iteration < 10; ++iteration) {
        size_t edges_before_tips = graph.edge_count();
        int tips     = remove_tips(graph, mean_read_length);
        graph.rebuild_adjacency();
        size_t edges_after_tips = graph.edge_count();
        size_t tip_edges_removed = edges_before_tips - edges_after_tips;

        size_t edges_before_internal = graph.edge_count();
        int internal_branches = remove_weak_internal_branches(graph, mean_read_length);
        graph.rebuild_adjacency();
        size_t edges_after_internal = graph.edge_count();
        size_t internal_edges_removed = edges_before_internal - edges_after_internal;

        size_t edges_before_low_wt = graph.edge_count();
        int low_wt   = prune_low_weight_edges(graph, options);
        graph.rebuild_adjacency();
        size_t edges_after_low_wt = graph.edge_count();
        size_t low_wt_edges_removed = edges_before_low_wt - edges_after_low_wt;

        spdlog::info(
            "Cleaning iteration {}: tip_nodes_removed={}, tip_edges_removed={}, "
            "internal_branch_nodes_removed={}, internal_branch_edges_removed={}, "
            "low_weight_edges_flagged={}, low_weight_edges_removed={}, "
            "bubble_popping_skipped=true, remaining_edges={}",
            iteration, tips, tip_edges_removed,
            internal_branches, internal_edges_removed,
            low_wt, low_wt_edges_removed,
            graph.edge_count());

        if (tips == 0 && internal_branches == 0 && low_wt == 0) break;
    }

    spdlog::info("Graph cleaning done: {} nodes, {} edges",
                 graph.node_count(), graph.edge_count());
}

} // namespace sharda
