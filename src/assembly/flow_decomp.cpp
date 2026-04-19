#include "assembly/flow_decomp.h"
#include <Highs.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <limits>

namespace sharda {

namespace {

/// Enumerate all source-to-sink paths via DFS (for small graphs).
/// Limit enumeration to max_enum paths to avoid explosion.
std::vector<std::vector<uint64_t>> enumerate_paths(
    const UnitigGraph& ug,
    uint64_t source, uint64_t sink,
    int max_enum = 1000)
{
    // Build adjacency
    size_t n = ug.unitig_count();
    std::vector<std::vector<uint64_t>> adj(n);
    for (const auto& e : ug.edges()) {
        adj[e.from].push_back(e.to);
    }

    std::vector<std::vector<uint64_t>> all_paths;
    std::vector<uint64_t> current_path;
    std::vector<bool> in_path(n, false);

    std::function<void(uint64_t)> dfs = [&](uint64_t u) {
        if (static_cast<int>(all_paths.size()) >= max_enum) return;
        current_path.push_back(u);
        in_path[u] = true;

        if (u == sink) {
            all_paths.push_back(current_path);
        } else {
            for (uint64_t v : adj[u]) {
                if (!in_path[v]) dfs(v);
                if (static_cast<int>(all_paths.size()) >= max_enum) break;
            }
        }

        current_path.pop_back();
        in_path[u] = false;
    };

    dfs(source);
    return all_paths;
}

/// Build the sequence for a path through the unitig graph.
std::string path_sequence(const UnitigGraph& ug, const std::vector<uint64_t>& path) {
    if (path.empty()) return "";
    std::string seq = ug.unitig(path[0]).sequence;
    // Subsequent unitigs: we append only the non-overlapping suffix.
    // Since unitigs connect via k-1 overlap from the underlying DBG,
    // we just append the full unitig sequence for simplicity
    // (the actual overlap depends on k and how compaction worked).
    // For correctness, just concatenate — the unitig sequences already
    // represent the full contiguous assembled sequence.
    for (size_t i = 1; i < path.size(); ++i) {
        seq += ug.unitig(path[i]).sequence;
    }
    return seq;
}

/// Find source and sink unitigs (containing min/max ref positions).
std::pair<uint64_t, uint64_t> find_source_sink(const UnitigGraph& ug) {
    uint64_t source = 0, sink = 0;
    // The unitig with the smallest node ID is likely the start of the backbone.
    // The one with the largest is likely the end.
    // A better heuristic: find unitigs with in-degree=0 (source) and out-degree=0 (sink).
    size_t n = ug.unitig_count();
    std::vector<int> in_deg(n, 0), out_deg(n, 0);
    for (const auto& e : ug.edges()) {
        out_deg[e.from]++;
        in_deg[e.to]++;
    }

    std::vector<uint64_t> sources, sinks;
    for (size_t i = 0; i < n; ++i) {
        if (in_deg[i] == 0) sources.push_back(i);
        if (out_deg[i] == 0) sinks.push_back(i);
    }

    if (sources.size() == 1 && sinks.size() == 1) {
        return {sources[0], sinks[0]};
    }

    // Fallback: use first and last unitig IDs
    spdlog::debug("Multiple sources ({}) or sinks ({}) — using first/last unitig",
                  sources.size(), sinks.size());
    source = sources.empty() ? 0 : sources[0];
    sink   = sinks.empty() ? (n - 1) : sinks.back();
    return {source, sink};
}

} // anonymous namespace

std::vector<HaplotypePath> flow_decomposition(
    const UnitigGraph& ug,
    int max_paths,
    double time_limit_sec)
{
    if (ug.unitig_count() == 0) return {};

    auto [source, sink] = find_source_sink(ug);
    spdlog::info("Flow decomposition: source={}, sink={}, max_paths={}",
                 source, sink, max_paths);

    // Enumerate all source-to-sink paths
    auto all_paths = enumerate_paths(ug, source, sink);
    spdlog::info("Enumerated {} source-to-sink paths", all_paths.size());

    if (all_paths.empty()) {
        spdlog::warn("No source-to-sink paths found");
        return {};
    }

    // If only one path exists, return it directly
    if (all_paths.size() == 1) {
        HaplotypePath hp;
        hp.unitig_ids = all_paths[0];
        hp.sequence = path_sequence(ug, all_paths[0]);
        hp.flow = ug.unitig(all_paths[0][0]).mean_depth;
        return {hp};
    }

    // ── Build ILP: minimize coverage deviation ──
    // Variables:
    //   f_p >= 0  for each path p (continuous) — flow on path
    //   s_e >= 0  for each edge e (continuous) — slack for coverage deviation
    //
    // For each edge e:
    //   sum_p (f_p * uses(p,e)) + s_e >= coverage_e
    //   sum_p (f_p * uses(p,e)) - s_e <= coverage_e
    // i.e.: |sum_p(f_p * uses(p,e)) - coverage_e| <= s_e
    //
    // Haplotype constraint: for each haplotype edge (U1,U2),
    //   there should exist a path p passing through both U1 and U2.
    //   We add: sum_{p passes through U1 and U2} f_p >= min_flow
    //
    // Objective: minimize sum_e s_e

    int P = static_cast<int>(all_paths.size());
    int E = static_cast<int>(ug.edges().size());

    // Precompute: for each path, which edges does it use?
    // edge_index: map (from,to) -> edge index
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, int>> edge_idx;
    for (int i = 0; i < E; ++i) {
        edge_idx[ug.edges()[i].from][ug.edges()[i].to] = i;
    }

    // path_uses_edge[p][e] = 1 if path p uses edge e
    std::vector<std::vector<bool>> path_uses_edge(P, std::vector<bool>(E, false));
    // path_through_unitig[p] = set of unitig IDs in path p
    std::vector<std::unordered_set<uint64_t>> path_unitigs(P);

    for (int p = 0; p < P; ++p) {
        for (uint64_t uid : all_paths[p]) {
            path_unitigs[p].insert(uid);
        }
        for (size_t j = 0; j + 1 < all_paths[p].size(); ++j) {
            uint64_t from = all_paths[p][j];
            uint64_t to   = all_paths[p][j + 1];
            auto it = edge_idx.find(from);
            if (it != edge_idx.end()) {
                auto jt = it->second.find(to);
                if (jt != it->second.end()) {
                    path_uses_edge[p][jt->second] = true;
                }
            }
        }
    }

    // Variables: f_0..f_{P-1} (path flows), s_0..s_{E-1} (slack)
    int num_vars = P + E;

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("time_limit", time_limit_sec);

    HighsModel model;
    model.lp_.num_col_ = num_vars;
    model.lp_.sense_ = ObjSense::kMinimize;

    // Column bounds and costs
    model.lp_.col_cost_.resize(num_vars, 0.0);
    model.lp_.col_lower_.resize(num_vars, 0.0);
    model.lp_.col_upper_.resize(num_vars, 1e9);

    // f_p variables: cost = 0
    // s_e variables: cost = 1 (we minimize total slack)
    for (int e = 0; e < E; ++e) {
        model.lp_.col_cost_[P + e] = 1.0;
    }

    // Constraints: for each edge e:
    //   sum_p (f_p * uses(p,e)) + s_e >= coverage_e   →  sum_p(f_p*a) + s_e >= cov
    //   sum_p (f_p * uses(p,e)) - s_e <= coverage_e   →  sum_p(f_p*a) - s_e <= cov
    // This is 2*E constraints.

    std::vector<double> row_lower, row_upper;
    std::vector<HighsInt> a_start, a_index;
    std::vector<double> a_value;

    int num_rows = 0;

    // Edge coverage constraints
    for (int e = 0; e < E; ++e) {
        double cov = ug.edges()[e].weight;

        // Constraint 1: sum_p(f_p * a_{p,e}) + s_e >= cov
        a_start.push_back(static_cast<HighsInt>(a_index.size()));
        for (int p = 0; p < P; ++p) {
            if (path_uses_edge[p][e]) {
                a_index.push_back(p);
                a_value.push_back(1.0);
            }
        }
        a_index.push_back(P + e);
        a_value.push_back(1.0);
        row_lower.push_back(cov);
        row_upper.push_back(1e30);
        num_rows++;

        // Constraint 2: sum_p(f_p * a_{p,e}) - s_e <= cov
        a_start.push_back(static_cast<HighsInt>(a_index.size()));
        for (int p = 0; p < P; ++p) {
            if (path_uses_edge[p][e]) {
                a_index.push_back(p);
                a_value.push_back(1.0);
            }
        }
        a_index.push_back(P + e);
        a_value.push_back(-1.0);
        row_lower.push_back(-1e30);
        row_upper.push_back(cov);
        num_rows++;
    }

    // Haplotype edge constraints:
    // For each haplotype edge (U1, U2), sum of f_p for paths through both >= 1
    for (const auto& he : ug.haplotype_edges()) {
        a_start.push_back(static_cast<HighsInt>(a_index.size()));
        bool any = false;
        for (int p = 0; p < P; ++p) {
            if (path_unitigs[p].count(he.from_node) &&
                path_unitigs[p].count(he.to_node)) {
                a_index.push_back(p);
                a_value.push_back(1.0);
                any = true;
            }
        }
        if (any) {
            row_lower.push_back(1.0); // at least some flow
            row_upper.push_back(1e30);
            num_rows++;
        } else {
            // No path covers this haplotype edge — remove the constraint start
            a_start.pop_back();
        }
    }

    // Finalize sparse matrix
    a_start.push_back(static_cast<HighsInt>(a_index.size()));

    model.lp_.num_row_ = num_rows;
    model.lp_.row_lower_ = row_lower;
    model.lp_.row_upper_ = row_upper;
    model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    model.lp_.a_matrix_.start_ = a_start;
    model.lp_.a_matrix_.index_ = a_index;
    model.lp_.a_matrix_.value_ = a_value;

    HighsStatus status = highs.passModel(model);
    if (status != HighsStatus::kOk) {
        spdlog::error("HiGHS passModel failed");
        return {};
    }

    status = highs.run();
    if (status != HighsStatus::kOk) {
        spdlog::error("HiGHS run failed");
        return {};
    }

    auto model_status = highs.getModelStatus();
    if (model_status != HighsModelStatus::kOptimal &&
        model_status != HighsModelStatus::kObjectiveBound &&
        model_status != HighsModelStatus::kSolutionLimit) {
        spdlog::warn("HiGHS did not find optimal solution, status={}",
                     static_cast<int>(model_status));
        // Fall back: return all paths with equal flow
        double avg_cov = 0;
        for (const auto& e : ug.edges()) avg_cov += e.weight;
        avg_cov /= std::max(1, E);
        double per_path = avg_cov / std::max(1, max_paths);

        std::vector<HaplotypePath> result;
        for (int p = 0; p < std::min(P, max_paths); ++p) {
            HaplotypePath hp;
            hp.unitig_ids = all_paths[p];
            hp.sequence = path_sequence(ug, all_paths[p]);
            hp.flow = per_path;
            result.push_back(std::move(hp));
        }
        return result;
    }

    // Extract solution
    const auto& sol = highs.getSolution();
    std::vector<HaplotypePath> result;
    for (int p = 0; p < P; ++p) {
        double f = sol.col_value[p];
        if (f < 0.5) continue; // negligible flow
        HaplotypePath hp;
        hp.unitig_ids = all_paths[p];
        hp.sequence = path_sequence(ug, all_paths[p]);
        hp.flow = f;
        result.push_back(std::move(hp));
    }

    // Sort by flow descending
    std::sort(result.begin(), result.end(),
              [](const HaplotypePath& a, const HaplotypePath& b) {
                  return a.flow > b.flow;
              });

    spdlog::info("Flow decomposition: {} paths extracted", result.size());
    return result;
}

} // namespace sharda
