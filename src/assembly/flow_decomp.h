#pragma once

#include "graph/unitig_graph.h"
#include <string>
#include <vector>
#include <utility>

namespace sharda {

struct HaplotypePath {
    std::vector<uint64_t> unitig_ids; // ordered unitig IDs in path
    std::string           sequence;
    double                flow;       // estimated coverage
};

/// Decompose the unitig graph into haplotype paths using ILP.
/// Returns empty vector on failure (e.g., infeasible).
/// max_paths: upper bound on number of paths (usually ploidy).
/// time_limit_sec: ILP solver time limit.
std::vector<HaplotypePath> flow_decomposition(
    const UnitigGraph& ug,
    int max_paths,
    double time_limit_sec = 60.0);

} // namespace sharda
