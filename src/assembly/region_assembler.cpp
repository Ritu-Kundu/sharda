#include "assembly/region_assembler.h"
#include "graph/dbg.h"
#include "graph/backbone.h"
#include "graph/unitig_graph.h"
#include "io/bam_reader.h"
#include "io/gfa_writer.h"
#include "assembly/read_adder.h"
#include "assembly/graph_cleaner.h"
#include "assembly/flow_decomp.h"
#include <spdlog/spdlog.h>
#include <filesystem>

namespace fs = std::filesystem;

namespace sharda {

RegionResult assemble_region(const RegionParams& params) {
    RegionResult result;
    result.region_name = params.region_name;

    try {
        // ── 1. Build backbone ───────────────────────────────────────
        spdlog::info("[{}] Building backbone (k={}, ref={} bp)",
                     params.region_name, params.k, params.ref_seq.size());
        DBG graph(params.k);
        build_backbone(graph, params.ref_seq, params.trs);

        // ── 2. Add reads (adjusting genomic → local coordinates) ────
        spdlog::info("[{}] Adding reads from: {}", params.region_name, params.bam_path);
        uint64_t read_pairs = 0;
        iterate_read_pairs(params.bam_path, [&](ReadPair&& pair) {
            // Convert genomic coordinates to local
            pair.read1.ref_start -= params.coord_offset;
            pair.read1.ref_end   -= params.coord_offset;
            pair.read2.ref_start -= params.coord_offset;
            pair.read2.ref_end   -= params.coord_offset;
            add_read_pair(pair, graph, params.trs);
            read_pairs++;
        });
        spdlog::info("[{}] Added {} read pairs, {} nodes, {} edges",
                     params.region_name, read_pairs,
                     graph.node_count(), graph.edge_count());

        // Debug GFA output
        if (!params.debug_dir.empty()) {
            fs::create_directories(params.debug_dir);
            write_gfa(params.debug_dir + "/raw.gfa", graph);
        }

        // ── 3. Clean graph ──────────────────────────────────────────
        spdlog::info("[{}] Cleaning graph", params.region_name);
        int mean_read_len = 150; // TODO: compute from actual reads
        clean_graph(graph, mean_read_len);

        if (!params.debug_dir.empty()) {
            write_gfa(params.debug_dir + "/clean.gfa", graph);
        }

        // ── 4. Build unitig graph ───────────────────────────────────
        spdlog::info("[{}] Building unitig graph", params.region_name);
        UnitigGraph ug;
        if (!ug.build(graph)) {
            result.error = "Unitig graph has cycles";
            spdlog::warn("[{}] {}", params.region_name, result.error);
            return result;
        }

        if (!params.debug_dir.empty()) {
            write_unitig_gfa(params.debug_dir + "/unitig.gfa", ug);
        }

        // ── 5. Flow decomposition ───────────────────────────────────
        spdlog::info("[{}] Flow decomposition (ploidy={})",
                     params.region_name, params.ploidy);
        auto paths = flow_decomposition(ug, params.ploidy);

        if (paths.empty()) {
            result.error = "No haplotype paths found";
            spdlog::warn("[{}] {}", params.region_name, result.error);
            return result;
        }

        // ── 6. Build output entries ─────────────────────────────────
        for (size_t i = 0; i < paths.size(); ++i) {
            std::string name = params.region_name + "_hap"
                             + std::to_string(i + 1)
                             + "_flow" + std::to_string(static_cast<int>(paths[i].flow));
            result.haplotypes.emplace_back(std::move(name), std::move(paths[i].sequence));
        }

        spdlog::info("[{}] Assembled {} haplotypes", params.region_name, paths.size());
        result.success = true;

    } catch (const std::exception& e) {
        result.error = e.what();
        spdlog::error("[{}] Failed: {}", params.region_name, result.error);
    }

    return result;
}

} // namespace sharda
