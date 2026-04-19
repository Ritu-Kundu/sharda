#include "assembly/anchor_chain.h"
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <climits>

namespace sharda {

std::vector<Anchor> find_and_chain_anchors(
    const std::vector<std::string>& read_kmers,
    const DBG& graph,
    int tr_id) {

    const auto& bb_nodes = graph.tr_nodes(tr_id);
    if (bb_nodes.empty()) return {};

    // Build kmer → list of backbone positions within this TR
    // Only keep kmers that appear exactly once (unique anchors).
    std::unordered_map<std::string, std::vector<uint64_t>> kmer_to_bb;
    for (uint64_t nid : bb_nodes) {
        const auto& nd = graph.node(nid);
        kmer_to_bb[nd.kmer].push_back(nid);
    }

    // Collect candidate anchors: (read_kmer_idx, backbone_node_id)
    std::vector<Anchor> candidates;
    for (int i = 0; i < static_cast<int>(read_kmers.size()); ++i) {
        auto it = kmer_to_bb.find(read_kmers[i]);
        if (it != kmer_to_bb.end() && it->second.size() == 1) {
            candidates.emplace_back(i, it->second[0]);
        }
    }

    if (candidates.empty()) return {};

    // Sort by read position (already in order, but be safe)
    std::sort(candidates.begin(), candidates.end());

    // DP chaining: find a subset of candidates that is co-linear in both
    // read position and backbone position, minimizing total gap discrepancy.
    // Simple O(n^2) DP — fine for read-length candidate lists.
    int n = static_cast<int>(candidates.size());
    std::vector<int64_t> dp(n, 0);   // best score ending at i
    std::vector<int>     prev(n, -1); // backtrack

    // Score: +1 per anchor, penalty for gap mismatch
    for (int i = 0; i < n; ++i) {
        dp[i] = 1; // just this anchor
        int ri = candidates[i].first;
        int32_t bi = graph.node(candidates[i].second).ref_pos;

        for (int j = 0; j < i; ++j) {
            int rj = candidates[j].first;
            int32_t bj = graph.node(candidates[j].second).ref_pos;

            // Must be co-linear: backbone pos must also increase
            if (bj >= bi) continue;

            int read_gap = ri - rj;
            int bb_gap   = bi - bj;
            int64_t penalty = std::abs(read_gap - bb_gap);
            int64_t score   = dp[j] + 1 - penalty;

            if (score > dp[i]) {
                dp[i]   = score;
                prev[i] = j;
            }
        }
    }

    // Backtrack from best endpoint
    int best_end = 0;
    for (int i = 1; i < n; ++i) {
        if (dp[i] > dp[best_end]) best_end = i;
    }

    std::vector<Anchor> chain;
    for (int i = best_end; i >= 0; i = prev[i]) {
        chain.push_back(candidates[i]);
        if (prev[i] < 0) break;
    }
    std::reverse(chain.begin(), chain.end());

    return chain;
}

} // namespace sharda
