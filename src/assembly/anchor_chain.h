#pragma once

#include "graph/dbg.h"
#include <string>
#include <vector>
#include <utility>

namespace sharda {

/// An anchor: (index of kmer in read, backbone node id in graph).
using Anchor = std::pair<int, uint64_t>;

/// Find unique-kmer anchors for a read within a specific TR's backbone nodes,
/// then chain them to minimize positional discrepancy.
/// Returns ordered anchors (read_kmer_idx, backbone_node_id).
std::vector<Anchor> find_and_chain_anchors(
    const std::vector<std::string>& read_kmers,
    const DBG& graph,
    int tr_id);

} // namespace sharda
