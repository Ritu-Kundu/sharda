#include "io/gfa_writer.h"
#include "io/debug_artifacts.h"
#include "graph/dbg.h"
#include "graph/unitig_graph.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <set>
#include <stdexcept>

namespace fs = std::filesystem;

namespace sharda {

namespace {

std::string json_escape(const std::string& value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (char ch : value) {
        switch (ch) {
        case '\\': escaped += "\\\\"; break;
        case '"': escaped += "\\\""; break;
        case '\n': escaped += "\\n"; break;
        case '\r': escaped += "\\r"; break;
        case '\t': escaped += "\\t"; break;
        default: escaped += ch; break;
        }
    }
    return escaped;
}

std::string stage_json_path(const DebugArtifactsConfig& config,
                            const std::string& stage_name) {
    return (fs::path(config.output_dir) / (stage_name + ".json")).string();
}

std::string stage_gfa_path(const DebugArtifactsConfig& config,
                           const std::string& stage_name) {
    return (fs::path(config.output_dir) / (stage_name + ".gfa")).string();
}

void write_viewer(const DebugArtifactsConfig& config) {
    if (!config.emit_html_view) {
        return;
    }

    std::ofstream out(fs::path(config.output_dir) / "viewer.html");
    if (!out) {
        throw std::runtime_error("Cannot open debug viewer in: " + config.output_dir);
    }

    out << R"HTML(<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Sharda Debug Viewer</title>
  <style>
    body { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; margin: 2rem; background: #f6f0e8; color: #1f1b17; }
    h1 { margin-bottom: 0.5rem; }
    p { max-width: 48rem; line-height: 1.5; }
    .hint { color: #5f564b; }
  </style>
</head>
<body>
  <h1>Sharda Debug Viewer</h1>
  <p class="hint">This initial viewer is intentionally static. Use the JSON files in this directory as the structured graph source.</p>
  <p class="hint">Expected files: raw.json, clean.json, unitig.json, manifest.json</p>
</body>
</html>
)HTML";
}

void write_manifest(const DebugArtifactsConfig& config,
                    const std::string& stage_name,
                    const std::string& graph_kind) {
    std::ofstream out(fs::path(config.output_dir) / "manifest.json");
    if (!out) {
        throw std::runtime_error("Cannot open debug manifest in: " + config.output_dir);
    }

    out << "{\n"
        << "  \"graph_format\": \"sharda-debug-v1\",\n"
        << "  \"latest_stage\": \"" << json_escape(stage_name) << "\",\n"
        << "  \"latest_graph_kind\": \"" << json_escape(graph_kind) << "\",\n"
        << "  \"stages\": [\n"
        << "    {\"name\": \"raw\", \"file\": \"raw.json\"},\n"
        << "    {\"name\": \"clean\", \"file\": \"clean.json\"},\n"
        << "    {\"name\": \"unitig\", \"file\": \"unitig.json\"}\n"
        << "  ]\n"
        << "}\n";
}

void write_read_trace_json_file(const std::string& path,
                                const std::vector<ReadTraceRecord>& traces) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open read trace JSON: " + path);
    }

    out << "{\n"
        << "  \"reads\": [\n";

    for (size_t trace_index = 0; trace_index < traces.size(); ++trace_index) {
        const auto& trace = traces[trace_index];
        out << "    {\"read_name\": \"" << json_escape(trace.read_name) << "\""
            << ", \"mate\": \"" << json_escape(trace.mate_label) << "\""
            << ", \"read_type\": \"" << json_escape(trace.read_type) << "\""
            << ", \"is_evidence\": " << (trace.is_evidence ? "true" : "false")
            << ", \"tr_id\": " << trace.tr_id
            << ", \"raw_nodes\": [";
        for (size_t node_index = 0; node_index < trace.raw_nodes.size(); ++node_index) {
            const auto& node = trace.raw_nodes[node_index];
            out << "{\"node_id\": " << node.node_id
                << ", \"sequence\": \"" << json_escape(node.sequence) << "\""
                << ", \"created\": " << (node.created ? "true" : "false")
                << ", \"is_backbone\": " << (node.is_backbone ? "true" : "false")
                << ", \"ref_pos\": " << node.ref_pos
                << ", \"tr_id\": " << node.tr_id
                << ", \"removed_after_clean\": " << (node.removed_after_clean ? "true" : "false");
            if (node.unitig_id == UINT64_MAX) {
                out << ", \"unitig_id\": null"
                    << ", \"unitig_sequence\": null";
            } else {
                out << ", \"unitig_id\": " << node.unitig_id
                    << ", \"unitig_sequence\": \"" << json_escape(node.unitig_sequence) << "\"";
            }
            out << "}";
            if (node_index + 1 != trace.raw_nodes.size()) {
                out << ", ";
            }
        }
        out << "]}";
        out << (trace_index + 1 == traces.size() ? "\n" : ",\n");
    }

    out << "  ]\n"
        << "}\n";
}

void write_locus_trace_json_file(const std::string& path,
                                 const std::vector<LocusTraceRecord>& traces) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open locus trace JSON: " + path);
    }

    out << "{\n"
        << "  \"loci\": [\n";

    for (size_t trace_index = 0; trace_index < traces.size(); ++trace_index) {
        const auto& trace = traces[trace_index];
        out << "    {\"local_start\": " << trace.local_start
            << ", \"length\": " << trace.length
            << ", \"global_start\": " << trace.global_start
            << ", \"reference_sequence\": \"" << json_escape(trace.reference_sequence) << "\""
            << ", \"raw_nodes\": [";
        for (size_t node_index = 0; node_index < trace.raw_nodes.size(); ++node_index) {
            const auto& node = trace.raw_nodes[node_index];
            out << "{\"node_id\": " << node.node_id
                << ", \"sequence\": \"" << json_escape(node.sequence) << "\""
                << ", \"is_backbone\": " << (node.is_backbone ? "true" : "false")
                << ", \"ref_pos\": " << node.ref_pos
                << ", \"tr_id\": " << node.tr_id
                << ", \"removed_after_clean\": " << (node.removed_after_clean ? "true" : "false");
            if (node.unitig_id == UINT64_MAX) {
                out << ", \"unitig_id\": null"
                    << ", \"unitig_sequence\": null";
            } else {
                out << ", \"unitig_id\": " << node.unitig_id
                    << ", \"unitig_sequence\": \"" << json_escape(node.unitig_sequence) << "\"";
            }
            out << "}";
            if (node_index + 1 != trace.raw_nodes.size()) {
                out << ", ";
            }
        }
        out << "]}";
        out << (trace_index + 1 == traces.size() ? "\n" : ",\n");
    }

    out << "  ]\n"
        << "}\n";
}

} // anonymous namespace

void write_gfa(const std::string& path, const DBG& graph) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open GFA: " + path);

    out << "H\tVN:Z:1.0\n";

    // S-lines (segments): each node
    for (const auto& node : graph.nodes()) {
        out << "S\t" << node.id << "\t" << node.kmer
            << "\tDP:f:" << node.depth
            << "\tBB:i:" << (node.is_backbone ? 1 : 0)
            << "\tRP:i:" << node.ref_pos
            << "\tTR:i:" << node.tr_id << '\n';
    }

    // L-lines (links): each edge
    for (const auto& edge : graph.edges()) {
        // Overlap = k-1
        int overlap = static_cast<int>(graph.k()) - 1;
        out << "L\t" << edge.from << "\t+\t" << edge.to << "\t+\t"
            << overlap << "M"
            << "\tRC:i:" << edge.weight << '\n';
    }
}

void write_dbg_json(const std::string& path, const DBG& graph) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open debug JSON: " + path);

    out << "{\n"
        << "  \"graph_kind\": \"dbg\",\n"
        << "  \"k\": " << graph.k() << ",\n"
        << "  \"nodes\": [\n";

    const auto& nodes = graph.nodes();
    for (size_t index = 0; index < nodes.size(); ++index) {
        const auto& node = nodes[index];
        out << "    {\"id\": " << node.id
            << ", \"sequence\": \"" << json_escape(node.kmer) << "\""
            << ", \"ref_pos\": " << node.ref_pos
            << ", \"tr_id\": " << node.tr_id
            << ", \"is_backbone\": " << (node.is_backbone ? "true" : "false")
            << ", \"depth\": " << node.depth
            << ", \"removed\": " << (graph.is_node_removed(node.id) ? "true" : "false")
            << "}";
        out << (index + 1 == nodes.size() ? "\n" : ",\n");
    }

    out << "  ],\n"
        << "  \"edges\": [\n";

    const auto& edges = graph.edges();
    for (size_t index = 0; index < edges.size(); ++index) {
        const auto& edge = edges[index];
        out << "    {\"from\": " << edge.from
            << ", \"to\": " << edge.to
            << ", \"weight\": " << edge.weight
            << "}";
        out << (index + 1 == edges.size() ? "\n" : ",\n");
    }

    out << "  ],\n"
        << "  \"haplotype_edges\": [\n";

    const auto& hap_edges = graph.haplotype_edges();
    for (size_t index = 0; index < hap_edges.size(); ++index) {
        const auto& edge = hap_edges[index];
        out << "    {\"from\": " << edge.from_node
            << ", \"to\": " << edge.to_node
            << ", \"weight\": " << edge.weight
            << "}";
        out << (index + 1 == hap_edges.size() ? "\n" : ",\n");
    }

    out << "  ]\n"
        << "}\n";
}

void write_unitig_gfa(const std::string& path, const UnitigGraph& ug) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open GFA: " + path);

    out << "H\tVN:Z:1.0\n";

    for (const auto& u : ug.unitigs()) {
        out << "S\t" << u.id << "\t" << u.sequence
            << "\tDP:f:" << u.mean_depth << '\n';
    }

    for (const auto& e : ug.edges()) {
        out << "L\t" << e.from << "\t+\t" << e.to << "\t+\t0M"
            << "\tRC:i:" << e.weight << '\n';
    }
}

void write_unitig_json(const std::string& path, const UnitigGraph& ug) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open debug JSON: " + path);

    out << "{\n"
        << "  \"graph_kind\": \"unitig\",\n"
        << "  \"unitigs\": [\n";

    const auto& unitigs = ug.unitigs();
    for (size_t index = 0; index < unitigs.size(); ++index) {
        const auto& unitig = unitigs[index];
        out << "    {\"id\": " << unitig.id
            << ", \"sequence\": \"" << json_escape(unitig.sequence) << "\""
            << ", \"mean_depth\": " << unitig.mean_depth
            << ", \"node_ids\": [";
        for (size_t node_index = 0; node_index < unitig.node_ids.size(); ++node_index) {
            out << unitig.node_ids[node_index];
            if (node_index + 1 != unitig.node_ids.size()) {
                out << ", ";
            }
        }
        out << "]}";
        out << (index + 1 == unitigs.size() ? "\n" : ",\n");
    }

    out << "  ],\n"
        << "  \"edges\": [\n";

    const auto& edges = ug.edges();
    for (size_t index = 0; index < edges.size(); ++index) {
        const auto& edge = edges[index];
        out << "    {\"from\": " << edge.from
            << ", \"to\": " << edge.to
            << ", \"weight\": " << edge.weight
            << "}";
        out << (index + 1 == edges.size() ? "\n" : ",\n");
    }

    out << "  ],\n"
        << "  \"haplotype_edges\": [\n";

    const auto& hap_edges = ug.haplotype_edges();
    for (size_t index = 0; index < hap_edges.size(); ++index) {
        const auto& edge = hap_edges[index];
        out << "    {\"from\": " << edge.from_node
            << ", \"to\": " << edge.to_node
            << ", \"weight\": " << edge.weight
            << "}";
        out << (index + 1 == hap_edges.size() ? "\n" : ",\n");
    }

    out << "  ]\n"
        << "}\n";
}

void ensure_debug_output_dir(const DebugArtifactsConfig& config) {
    if (!config.should_write()) {
        return;
    }

    fs::create_directories(config.output_dir);
    write_viewer(config);
}

void write_dbg_debug_artifacts(const DebugArtifactsConfig& config,
                               const std::string& stage_name,
                               const DBG& graph) {
    if (!config.should_write()) {
        return;
    }

    ensure_debug_output_dir(config);
    if (config.emit_gfa) {
        write_gfa(stage_gfa_path(config, stage_name), graph);
    }
    if (config.emit_json) {
        write_dbg_json(stage_json_path(config, stage_name), graph);
    }
    write_manifest(config, stage_name, "dbg");
}

void write_unitig_debug_artifacts(const DebugArtifactsConfig& config,
                                  const std::string& stage_name,
                                  const UnitigGraph& graph) {
    if (!config.should_write()) {
        return;
    }

    ensure_debug_output_dir(config);
    if (config.emit_gfa) {
        write_unitig_gfa(stage_gfa_path(config, stage_name), graph);
    }
    if (config.emit_json) {
        write_unitig_json(stage_json_path(config, stage_name), graph);
    }
    write_manifest(config, stage_name, "unitig");
}

std::vector<LocusTraceRecord> collect_locus_traces(
    const std::vector<LocusTraceRequest>& requests,
    const std::string& reference_sequence,
    int32_t coord_offset,
    const DBG& graph) {
    std::vector<LocusTraceRecord> traces;
    traces.reserve(requests.size());

    for (const auto& request : requests) {
        int32_t start = std::max<int32_t>(0, request.start);
        int32_t length = std::max<int32_t>(0, request.length);
        int32_t end = std::min<int32_t>(static_cast<int32_t>(reference_sequence.size()), start + length);

        LocusTraceRecord trace;
        trace.local_start = start;
        trace.length = end - start;
        trace.global_start = coord_offset + start;
        trace.reference_sequence = reference_sequence.substr(start, trace.length);

        std::set<uint64_t> node_ids;
        for (const auto& node : graph.nodes()) {
            if (node.ref_pos >= start && node.ref_pos < end) {
                node_ids.insert(node.id);
            }
        }

        std::set<uint64_t> expanded_ids = node_ids;
        for (uint64_t node_id : node_ids) {
            for (uint64_t edge_index : graph.out_edges(node_id)) {
                expanded_ids.insert(graph.edges()[edge_index].to);
            }
            for (uint64_t edge_index : graph.in_edges(node_id)) {
                expanded_ids.insert(graph.edges()[edge_index].from);
            }
        }

        for (uint64_t node_id : expanded_ids) {
            const auto& node = graph.node(node_id);
            trace.raw_nodes.push_back({
                node.id,
                node.kmer,
                node.is_backbone,
                node.ref_pos,
                node.tr_id,
                false,
                UINT64_MAX,
                ""
            });
        }

        std::sort(trace.raw_nodes.begin(), trace.raw_nodes.end(),
                  [](const LocusTraceNode& left, const LocusTraceNode& right) {
                      if (left.ref_pos != right.ref_pos) {
                          return left.ref_pos < right.ref_pos;
                      }
                      return left.node_id < right.node_id;
                  });

        traces.push_back(std::move(trace));
    }

    return traces;
}

void finalize_locus_traces(std::vector<LocusTraceRecord>& traces,
                           const DBG& graph,
                           const UnitigGraph& unitig_graph) {
    for (auto& trace : traces) {
        for (auto& node : trace.raw_nodes) {
            node.removed_after_clean = graph.is_node_removed(node.node_id);
            if (node.removed_after_clean) {
                continue;
            }

            node.unitig_id = unitig_graph.node_to_unitig(node.node_id);
            if (node.unitig_id != UINT64_MAX) {
                node.unitig_sequence = unitig_graph.unitig(node.unitig_id).sequence;
            }
        }
    }
}

void write_locus_trace_artifacts(const DebugArtifactsConfig& config,
                                 const std::vector<LocusTraceRecord>& traces) {
    if (!config.should_trace_loci() || traces.empty()) {
        return;
    }

    ensure_debug_output_dir(config);
    write_locus_trace_json_file((fs::path(config.output_dir) / "locus_traces.json").string(),
                                traces);
}

void finalize_read_traces(std::vector<ReadTraceRecord>& traces,
                          const DBG& graph,
                          const UnitigGraph& unitig_graph) {
    for (auto& trace : traces) {
        for (auto& node : trace.raw_nodes) {
            node.removed_after_clean = graph.is_node_removed(node.node_id);
            if (node.removed_after_clean) {
                node.unitig_id = UINT64_MAX;
                node.unitig_sequence.clear();
                continue;
            }

            node.unitig_id = unitig_graph.node_to_unitig(node.node_id);
            if (node.unitig_id != UINT64_MAX) {
                node.unitig_sequence = unitig_graph.unitig(node.unitig_id).sequence;
            }
        }
    }
}

void write_read_trace_artifacts(const DebugArtifactsConfig& config,
                                const std::vector<ReadTraceRecord>& traces) {
    if (!config.should_trace_reads() || traces.empty()) {
        return;
    }

    ensure_debug_output_dir(config);
    write_read_trace_json_file((fs::path(config.output_dir) / "read_traces.json").string(),
                               traces);
}

} // namespace sharda
