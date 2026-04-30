#include "graph/backbone.h"
#include "util/kmer.h"
#include <spdlog/spdlog.h>
#include <algorithm>
#include <chrono>

namespace sharda {

void build_backbone(DBG& graph,
                    const std::string& ref_seq,
                    const std::vector<TandemRepeat>& trs) {
    const auto stage_start = std::chrono::steady_clock::now();
    const auto elapsed_ms = [&]() {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::steady_clock::now() - stage_start)
            .count();
    };

    int k = graph.k();
    if (static_cast<int>(ref_seq.size()) < k) {
        spdlog::warn("Reference shorter than k={}, no backbone built", k);
        return;
    }

    spdlog::info("Backbone stage: extracting {}-mers from {} bp reference", k, ref_seq.size());
    auto kmers = extract_kmers(ref_seq, k);
    spdlog::info("Backbone stage: extracted {} kmers ({} ms)", kmers.size(), elapsed_ms());

    // Pre-build sorted TR intervals for quick lookup
    // For each ref_pos, determine which TR (if any) it falls within.
    // A kmer at ref_pos covers [ref_pos, ref_pos+k). We say it belongs to
    // a TR if its start position (ref_pos) is within the TR interval.
    auto tr_id_for_pos = [&](int32_t pos) -> int {
        for (const auto& tr : trs) {
            if (pos >= tr.start && pos < tr.end)
                return tr.id;
        }
        return -1;
    };

    uint64_t prev_id = UINT64_MAX;
    size_t next_progress = 250;
    for (size_t i = 0; i < kmers.size(); ++i) {
        int32_t ref_pos = static_cast<int32_t>(i);
        int tr_id = tr_id_for_pos(ref_pos);
        uint64_t nid = graph.add_backbone_node(kmers[i], ref_pos, tr_id);

        if (prev_id != UINT64_MAX) {
            graph.add_edge(prev_id, nid);
        }
        prev_id = nid;

        if (i + 1 == next_progress || i + 1 == kmers.size()) {
            spdlog::info(
                "Backbone stage: inserted {}/{} kmers (nodes={}, edges={}, {} ms)",
                i + 1, kmers.size(), graph.node_count(), graph.edge_count(), elapsed_ms());
            next_progress += 250;
        }
    }

    spdlog::info("Backbone built: {} nodes, {} edges, k={}",
                 graph.node_count(), graph.edge_count(), k);
}

} // namespace sharda
