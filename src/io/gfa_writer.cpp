#include "io/gfa_writer.h"
#include "io/debug_artifacts.h"
#include "graph/dbg.h"
#include "graph/unitig_graph.h"
#include <filesystem>
#include <fstream>
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

} // namespace sharda
