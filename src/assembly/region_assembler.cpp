#include "assembly/region_assembler.h"
#include "graph/dbg.h"
#include "graph/backbone.h"
#include "graph/unitig_graph.h"
#include "assembly/sv_caller.h"
#include "io/bam_reader.h"
#include "io/debug_artifacts.h"
#include "io/gfa_writer.h"
#include "assembly/read_adder.h"
#include "assembly/graph_cleaner.h"
#include "assembly/flow_decomp.h"
#include <spdlog/spdlog.h>
#include <algorithm>
#include <filesystem>
#include <functional>
#include <limits>
#include <unordered_map>
#include <unordered_set>

namespace fs = std::filesystem;

namespace sharda {

namespace {

std::string infer_chrom_name(const std::string& region_name) {
    const size_t separator = region_name.find(':');
    if (separator == std::string::npos) {
        return region_name;
    }
    return region_name.substr(0, separator);
}

std::vector<std::vector<uint64_t>> enumerate_unitig_paths(
    const UnitigGraph& graph,
    uint64_t source,
    uint64_t sink,
    size_t max_paths) {
    std::vector<std::vector<uint64_t>> adjacency(graph.unitig_count());
    for (const auto& edge : graph.edges()) {
        adjacency[edge.from].push_back(edge.to);
    }

    std::vector<std::vector<uint64_t>> paths;
    std::vector<uint64_t> current_path;
    std::vector<bool> in_path(graph.unitig_count(), false);

    std::function<void(uint64_t)> dfs = [&](uint64_t unitig_id) {
        if (paths.size() >= max_paths) {
            return;
        }

        current_path.push_back(unitig_id);
        in_path[unitig_id] = true;

        if (unitig_id == sink) {
            paths.push_back(current_path);
        } else {
            for (uint64_t next : adjacency[unitig_id]) {
                if (!in_path[next]) {
                    dfs(next);
                    if (paths.size() >= max_paths) {
                        break;
                    }
                }
            }
        }

        current_path.pop_back();
        in_path[unitig_id] = false;
    };

    dfs(source);
    return paths;
}

struct AlternateIntervalPath {
    uint64_t sink_id;
    std::vector<uint64_t> path;
};

std::vector<AlternateIntervalPath> enumerate_alternate_interval_paths(
    const UnitigGraph& graph,
    const std::vector<std::vector<uint64_t>>& adjacency,
    uint64_t source,
    size_t max_paths) {
    std::vector<AlternateIntervalPath> paths;
    std::vector<uint64_t> current_path{source};
    std::vector<bool> in_path(graph.unitig_count(), false);
    in_path[source] = true;

    std::function<void(uint64_t)> dfs = [&](uint64_t unitig_id) {
        if (paths.size() >= max_paths || in_path[unitig_id]) {
            return;
        }

        current_path.push_back(unitig_id);
        in_path[unitig_id] = true;

        const auto& unitig = graph.unitig(unitig_id);
        if (unitig.support_class() == UnitigSupportClass::BackboneOnly) {
            paths.push_back({unitig_id, current_path});
        } else {
            for (uint64_t next : adjacency[unitig_id]) {
                dfs(next);
                if (paths.size() >= max_paths) {
                    break;
                }
            }
        }

        in_path[unitig_id] = false;
        current_path.pop_back();
    };

    for (uint64_t next : adjacency[source]) {
        if (graph.unitig(next).support_class() == UnitigSupportClass::BackboneOnly) {
            continue;
        }
        dfs(next);
        if (paths.size() >= max_paths) {
            break;
        }
    }

    return paths;
}

bool is_backbone_only_path(const UnitigGraph& graph,
                          const std::vector<uint64_t>& path) {
    return std::all_of(path.begin(), path.end(), [&](uint64_t unitig_id) {
        return graph.unitig(unitig_id).support_class() == UnitigSupportClass::BackboneOnly;
    });
}

bool path_has_non_backbone_interval(const UnitigGraph& graph,
                                    const std::vector<uint64_t>& path) {
    if (path.size() <= 2) {
        return false;
    }
    for (size_t index = 1; index + 1 < path.size(); ++index) {
        if (graph.unitig(path[index]).support_class() != UnitigSupportClass::BackboneOnly) {
            return true;
        }
    }
    return false;
}

std::string reconstruct_path_sequence(const UnitigGraph& graph,
                                      const std::vector<uint64_t>& path) {
    if (path.empty()) {
        return "";
    }

    std::string sequence = graph.unitig(path.front()).sequence;
    const size_t overlap = graph.k() > 0 ? static_cast<size_t>(graph.k() - 1) : 0;
    for (size_t index = 1; index < path.size(); ++index) {
        const auto& next_sequence = graph.unitig(path[index]).sequence;
        if (next_sequence.size() <= overlap) {
            continue;
        }
        sequence += next_sequence.substr(overlap);
    }

    return sequence;
}

bool build_indel_call(const std::string& chrom,
                      int32_t coord_offset,
                      uint64_t source_unitig_id,
                      uint64_t sink_unitig_id,
                      int32_t source_ref_pos,
                      int32_t sink_ref_pos,
                      int32_t reference_start,
                      const std::string& reference_sequence,
                      const std::string& alternate_sequence,
                      const std::string& call_id,
                      StructuralVariantCall& call) {
    if (reference_sequence.empty() || alternate_sequence.empty()
        || reference_sequence == alternate_sequence) {
        return false;
    }

    size_t prefix = 0;
    while (prefix < reference_sequence.size()
           && prefix < alternate_sequence.size()
           && reference_sequence[prefix] == alternate_sequence[prefix]) {
        ++prefix;
    }

    size_t ref_suffix = reference_sequence.size();
    size_t alt_suffix = alternate_sequence.size();
    while (ref_suffix > prefix && alt_suffix > prefix
           && reference_sequence[ref_suffix - 1] == alternate_sequence[alt_suffix - 1]) {
        --ref_suffix;
        --alt_suffix;
    }

    if (prefix == 0) {
        return false;
    }

    const std::string ref_core = reference_sequence.substr(prefix, ref_suffix - prefix);
    const std::string alt_core = alternate_sequence.substr(prefix, alt_suffix - prefix);
    if (!ref_core.empty() && !alt_core.empty()) {
        return false;
    }

    const size_t anchor_index = prefix - 1;
    const std::string ref_allele = reference_sequence.substr(anchor_index, ref_suffix - anchor_index);
    const std::string alt_allele = alternate_sequence.substr(anchor_index, alt_suffix - anchor_index);
    if (ref_allele.empty() || alt_allele.empty()) {
        return false;
    }

    call.chrom = chrom;
    call.pos = coord_offset + reference_start + static_cast<int32_t>(anchor_index);
    call.end = call.pos + static_cast<int32_t>(ref_allele.size());
    call.id = call_id;
    call.ref = ref_allele;
    call.alt = alt_allele;
    call.sv_type = alt_allele.size() > ref_allele.size() ? "INS" : "DEL";
    call.sv_len = static_cast<int32_t>(alt_allele.size()) - static_cast<int32_t>(ref_allele.size());
    call.info_fields.push_back("SRC_UID=" + std::to_string(source_unitig_id));
    call.info_fields.push_back("SNK_UID=" + std::to_string(sink_unitig_id));
    call.info_fields.push_back("SRC_REF_POS=" + std::to_string(coord_offset + source_ref_pos + 1));
    call.info_fields.push_back("SNK_REF_POS=" + std::to_string(coord_offset + sink_ref_pos + 1));
    return true;
}

double score_alternate_path(const UnitigGraph& graph,
                            const std::vector<uint64_t>& path) {
    double min_non_backbone_depth = std::numeric_limits<double>::infinity();
    for (size_t index = 1; index + 1 < path.size(); ++index) {
        const auto& unitig = graph.unitig(path[index]);
        if (unitig.support_class() == UnitigSupportClass::BackboneOnly) {
            continue;
        }
        min_non_backbone_depth = std::min(min_non_backbone_depth, unitig.mean_depth);
    }

    if (!std::isfinite(min_non_backbone_depth)) {
        return 0.0;
    }
    return min_non_backbone_depth;
}

std::string collapse_key(const StructuralVariantCall& call) {
    return call.chrom + "\n"
        + std::to_string(call.pos) + "\n"
        + std::to_string(call.end) + "\n"
        + call.ref + "\n"
        + call.alt + "\n"
        + call.sv_type + "\n"
        + std::to_string(call.sv_len);
}

int32_t extract_info_int(const StructuralVariantCall& call,
                         const std::string& prefix) {
    for (const auto& field : call.info_fields) {
        if (field.rfind(prefix, 0) == 0) {
            return std::stoi(field.substr(prefix.size()));
        }
    }
    return -1;
}

bool call_has_better_rank(const StructuralVariantCall& candidate,
                          const StructuralVariantCall& incumbent) {
    if (candidate.support_score != incumbent.support_score) {
        return candidate.support_score > incumbent.support_score;
    }

    const int32_t candidate_span = extract_info_int(candidate, "SNK_REF_POS=")
        - extract_info_int(candidate, "SRC_REF_POS=");
    const int32_t incumbent_span = extract_info_int(incumbent, "SNK_REF_POS=")
        - extract_info_int(incumbent, "SRC_REF_POS=");
    if (candidate_span != incumbent_span) {
        return candidate_span < incumbent_span;
    }

    if (candidate.chrom != incumbent.chrom) {
        return candidate.chrom < incumbent.chrom;
    }
    if (candidate.pos != incumbent.pos) {
        return candidate.pos < incumbent.pos;
    }
    if (candidate.end != incumbent.end) {
        return candidate.end < incumbent.end;
    }
    return candidate.alt < incumbent.alt;
}

bool intervals_overlap_or_touch(const StructuralVariantCall& lhs,
                                const StructuralVariantCall& rhs) {
    const int32_t lhs_start = lhs.pos;
    const int32_t lhs_end = std::max(lhs.pos, lhs.end);
    const int32_t rhs_start = rhs.pos;
    const int32_t rhs_end = std::max(rhs.pos, rhs.end);
    return lhs_start <= rhs_end && rhs_start <= lhs_end;
}

bool should_merge_overlapping_calls(const StructuralVariantCall& lhs,
                                    const StructuralVariantCall& rhs) {
    if (lhs.chrom != rhs.chrom || lhs.sv_type != rhs.sv_type || lhs.sv_len != rhs.sv_len) {
        return false;
    }
    if (lhs.sv_type != "DEL") {
        return false;
    }
    if (lhs.alt != rhs.alt) {
        return false;
    }
    if (!intervals_overlap_or_touch(lhs, rhs)) {
        return false;
    }

    const int32_t lhs_source_ref = extract_info_int(lhs, "SRC_REF_POS=");
    const int32_t lhs_sink_ref = extract_info_int(lhs, "SNK_REF_POS=");
    const int32_t rhs_source_ref = extract_info_int(rhs, "SRC_REF_POS=");
    const int32_t rhs_sink_ref = extract_info_int(rhs, "SNK_REF_POS=");

    return lhs_source_ref == rhs_source_ref || lhs_sink_ref == rhs_sink_ref;
}

} // namespace

std::vector<StructuralVariantCall> collapse_structural_variant_calls(
    const std::vector<StructuralVariantCall>& calls) {
    std::vector<StructuralVariantCall> collapsed;
    std::unordered_map<std::string, size_t> best_call_by_key;

    for (const auto& call : calls) {
        const std::string key = collapse_key(call);
        auto it = best_call_by_key.find(key);
        if (it == best_call_by_key.end()) {
            best_call_by_key.emplace(key, collapsed.size());
            collapsed.push_back(call);
            continue;
        }

        if (call_has_better_rank(call, collapsed[it->second])) {
            collapsed[it->second] = call;
        }
    }

    bool merged_overlap = true;
    while (merged_overlap) {
        merged_overlap = false;
        for (size_t i = 0; i < collapsed.size() && !merged_overlap; ++i) {
            for (size_t j = i + 1; j < collapsed.size(); ++j) {
                if (!should_merge_overlapping_calls(collapsed[i], collapsed[j])) {
                    continue;
                }

                if (call_has_better_rank(collapsed[j], collapsed[i])) {
                    collapsed[i] = collapsed[j];
                }
                collapsed.erase(collapsed.begin() + static_cast<std::ptrdiff_t>(j));
                merged_overlap = true;
                break;
            }
        }
    }

    std::sort(collapsed.begin(), collapsed.end(),
              [](const StructuralVariantCall& lhs, const StructuralVariantCall& rhs) {
                  if (lhs.chrom != rhs.chrom) {
                      return lhs.chrom < rhs.chrom;
                  }
                  if (lhs.pos != rhs.pos) {
                      return lhs.pos < rhs.pos;
                  }
                  if (lhs.end != rhs.end) {
                      return lhs.end < rhs.end;
                  }
                  return lhs.alt < rhs.alt;
              });

    for (size_t index = 0; index < collapsed.size(); ++index) {
        collapsed[index].id = "sv" + std::to_string(index + 1);
    }

    return collapsed;
}

std::vector<StructuralVariantCall> call_structural_variants(
    const UnitigGraph& graph,
    const std::string& chrom,
    int32_t coord_offset,
    size_t max_paths_per_interval) {
    std::vector<StructuralVariantCall> calls;
    if (graph.unitig_count() == 0) {
        return calls;
    }

    std::vector<int> in_degree(graph.unitig_count(), 0);
    std::vector<int> out_degree(graph.unitig_count(), 0);
    std::vector<std::vector<uint64_t>> adjacency(graph.unitig_count());
    for (const auto& edge : graph.edges()) {
        adjacency[edge.from].push_back(edge.to);
        out_degree[edge.from]++;
        in_degree[edge.to]++;
    }

    size_t call_index = 0;

    for (size_t source_id = 0; source_id < graph.unitig_count(); ++source_id) {
        const auto& source = graph.unitig(source_id);
        if (source.support_class() != UnitigSupportClass::BackboneOnly || out_degree[source_id] <= 1) {
            continue;
        }

        auto alternate_interval_paths = enumerate_alternate_interval_paths(
            graph,
            adjacency,
            source_id,
            max_paths_per_interval);

        for (const auto& alternate_interval : alternate_interval_paths) {
            const uint64_t sink_id = alternate_interval.sink_id;
            if (sink_id == source_id) {
                continue;
            }

            const auto& sink = graph.unitig(sink_id);
            if (sink.support_class() != UnitigSupportClass::BackboneOnly || in_degree[sink_id] <= 1) {
                continue;
            }
            if (sink.ref_pos <= source.ref_pos) {
                continue;
            }

            auto paths = enumerate_unitig_paths(graph, source_id, sink_id, max_paths_per_interval);
            if (paths.size() < 2) {
                continue;
            }

            std::vector<std::vector<uint64_t>> reference_paths;
            for (const auto& path : paths) {
                if (is_backbone_only_path(graph, path)) {
                    reference_paths.push_back(path);
                }
            }

            if (reference_paths.size() != 1 || !path_has_non_backbone_interval(graph, alternate_interval.path)) {
                continue;
            }

            const std::string reference_sequence = reconstruct_path_sequence(graph, reference_paths.front());
            StructuralVariantCall call;
            const std::string alternate_sequence = reconstruct_path_sequence(graph, alternate_interval.path);
            if (!build_indel_call(chrom,
                                  coord_offset,
                                  source_id,
                                  sink_id,
                                  source.ref_pos,
                                  sink.ref_pos,
                                  source.ref_pos,
                                  reference_sequence,
                                  alternate_sequence,
                                  "sv" + std::to_string(++call_index),
                                  call)) {
                continue;
            }
            call.support_score = score_alternate_path(graph, alternate_interval.path);
            calls.push_back(std::move(call));
        }
    }

    return collapse_structural_variant_calls(calls);
}

RegionResult assemble_region(const RegionParams& params) {
    RegionResult result;
    result.region_name = params.region_name;
    std::vector<ReadTraceRecord> read_traces;
    ReadTraceSink trace_sink{&params.debug_artifacts.traced_reads, &read_traces};
    std::vector<LocusTraceRecord> locus_traces;

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
            add_read_pair(pair, graph, params.trs, trace_sink);
            read_pairs++;
        });
        spdlog::info("[{}] Added {} read pairs, {} nodes, {} edges",
                     params.region_name, read_pairs,
                     graph.node_count(), graph.edge_count());

        if (params.debug_artifacts.should_trace_loci()) {
            locus_traces = collect_locus_traces(params.debug_artifacts.traced_loci,
                                                params.ref_seq,
                                                params.coord_offset,
                                                graph);
        }

        if (params.debug_artifacts.should_write()) {
            write_dbg_debug_artifacts(params.debug_artifacts, "raw", graph);
        } else if (!params.debug_dir.empty()) {
            fs::create_directories(params.debug_dir);
            write_gfa(params.debug_dir + "/raw.gfa", graph);
        }

        // ── 3. Clean graph ──────────────────────────────────────────
        spdlog::info("[{}] Cleaning graph", params.region_name);
        int mean_read_len = 150; // TODO: compute from actual reads
        GraphCleaningOptions cleaning_options;
        cleaning_options.preserve_backbone_edges = execution_mode_runs_sv(params.mode);
        clean_graph(graph, mean_read_len, cleaning_options);

        if (params.debug_artifacts.should_write()) {
            write_dbg_debug_artifacts(params.debug_artifacts, "clean", graph);
        } else if (!params.debug_dir.empty()) {
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

        if (params.debug_artifacts.should_trace_reads()) {
            finalize_read_traces(read_traces, graph, ug);
            write_read_trace_artifacts(params.debug_artifacts, read_traces);
        }
        if (params.debug_artifacts.should_trace_loci()) {
            finalize_locus_traces(locus_traces, graph, ug);
            write_locus_trace_artifacts(params.debug_artifacts, locus_traces);
        }

        if (params.debug_artifacts.should_write()) {
            write_unitig_debug_artifacts(params.debug_artifacts, "unitig", ug);
            if (execution_mode_runs_sv(params.mode)) {
                write_sv_unitig_gfa((fs::path(params.debug_artifacts.output_dir) / "unitig.sv.gfa").string(), ug);
                write_sv_unitig_json((fs::path(params.debug_artifacts.output_dir) / "unitig.sv.json").string(), ug);
            }
        } else if (!params.debug_dir.empty()) {
            write_unitig_gfa(params.debug_dir + "/unitig.gfa", ug);
            if (execution_mode_runs_sv(params.mode)) {
                write_sv_unitig_gfa(params.debug_dir + "/unitig.sv.gfa", ug);
                write_sv_unitig_json(params.debug_dir + "/unitig.sv.json", ug);
            }
        }

        if (params.stop_after_unitig_graph) {
            spdlog::info("[{}] Stopping after unitig graph construction", params.region_name);
            result.success = true;
            return result;
        }

        if (execution_mode_runs_sv(params.mode)) {
            result.sv_calls = call_structural_variants(
                ug,
                infer_chrom_name(params.region_name),
                params.coord_offset);
        }

        if (!execution_mode_runs_haplotype(params.mode)) {
            spdlog::info("[{}] Skipping haplotype flow decomposition in SV-only mode",
                         params.region_name);
            result.success = true;
            return result;
        }

        // ── 5. Flow decomposition ───────────────────────────────────
        spdlog::info("[{}] Flow decomposition (ploidy={})",
                     params.region_name, params.ploidy);
        FlowBoundaryAnchors anchors;
        anchors.start_node_id = graph.backbone_node_at(0);
        if (static_cast<int>(params.ref_seq.size()) >= params.k) {
            anchors.end_node_id = graph.backbone_node_at(
                static_cast<int32_t>(params.ref_seq.size()) - params.k);
        }
        auto paths = flow_decomposition(ug, params.ploidy, anchors);

        if (params.debug_artifacts.should_write()) {
            write_flow_path_artifacts(params.debug_artifacts, paths);
        }

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
