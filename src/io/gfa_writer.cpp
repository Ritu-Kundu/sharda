#include "io/gfa_writer.h"
#include "io/debug_artifacts.h"
#include "io/vcf_writer.h"
#include "graph/dbg.h"
#include "graph/unitig_graph.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
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

template <typename Container>
void write_ref_positions_json(std::ostream& out, const Container& positions) {
    out << '[';
    bool first = true;
    for (int32_t pos : positions) {
        if (!first) {
            out << ", ";
        }
        out << pos;
        first = false;
    }
    out << ']';
}

std::string ref_positions_tag(const std::set<int32_t>& positions) {
    if (positions.empty()) {
        return "";
    }

    std::string tag;
    bool first = true;
    for (int32_t pos : positions) {
        if (!first) {
            tag += ',';
        }
        tag += std::to_string(pos);
        first = false;
    }
    return tag;
}

std::string stage_json_path(const DebugArtifactsConfig& config,
                            const std::string& stage_name) {
    return (fs::path(config.output_dir) / (stage_name + ".json")).string();
}

std::string stage_gfa_path(const DebugArtifactsConfig& config,
                           const std::string& stage_name) {
    return (fs::path(config.output_dir) / (stage_name + ".gfa")).string();
}

std::string flow_paths_json_path(const DebugArtifactsConfig& config) {
    return (fs::path(config.output_dir) / "flow_paths.json").string();
}

std::string unitig_support_color(const Unitig& unitig) {
    switch (unitig.support_class()) {
    case UnitigSupportClass::BackboneOnly:
        return "#3B7A57";
    case UnitigSupportClass::Mixed:
        return "#C97B2B";
    case UnitigSupportClass::ReadOnly:
        return "#2D6A9F";
    }
    return "#2D6A9F";
}

std::string escape_vcf_field(const std::string& value) {
    return value.empty() ? "." : value;
}

std::string format_vcf_float(double value) {
    std::ostringstream out;
    out << std::setprecision(6) << std::defaultfloat << value;
    return out.str();
}

std::string build_vcf_info_field(const StructuralVariantCall& call) {
    std::vector<std::string> info_fields = call.info_fields;

    auto append_if_missing = [&info_fields](const std::string& prefix,
                                            const std::string& value) {
        for (const auto& field : info_fields) {
            if (field.rfind(prefix, 0) == 0) {
                return;
            }
        }
        info_fields.push_back(prefix + value);
    };

    append_if_missing("SVTYPE=", call.sv_type.empty() ? "UNK" : call.sv_type);
    append_if_missing("END=", std::to_string(call.end));
    append_if_missing("SVLEN=", std::to_string(call.sv_len));
    append_if_missing("SUPPORT=", format_vcf_float(call.support_score));

    if (info_fields.empty()) {
        return ".";
    }

    std::ostringstream out;
    for (size_t index = 0; index < info_fields.size(); ++index) {
        if (index > 0) {
            out << ';';
        }
        out << info_fields[index];
    }
    return out.str();
}

struct StableNodeOrder {
    std::vector<uint64_t> ordered_internal_ids;
    std::vector<uint64_t> stable_id_by_internal;
};

bool ref_positions_less(const std::set<int32_t>& lhs,
                        const std::set<int32_t>& rhs) {
    return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

StableNodeOrder build_stable_node_order(const DBG& graph) {
    StableNodeOrder order;
    order.ordered_internal_ids.resize(graph.node_count());
    std::iota(order.ordered_internal_ids.begin(), order.ordered_internal_ids.end(), 0);

    std::stable_sort(order.ordered_internal_ids.begin(), order.ordered_internal_ids.end(),
                     [&graph](uint64_t lhs_id, uint64_t rhs_id) {
        const auto& lhs = graph.node(lhs_id);
        const auto& rhs = graph.node(rhs_id);

        if (lhs.is_backbone != rhs.is_backbone) {
            return lhs.is_backbone > rhs.is_backbone;
        }
        if (lhs.ref_pos != rhs.ref_pos) {
            return lhs.ref_pos < rhs.ref_pos;
        }
        if (lhs.kmer != rhs.kmer) {
            return lhs.kmer < rhs.kmer;
        }
        if (ref_positions_less(lhs.ref_positions, rhs.ref_positions)) {
            return true;
        }
        if (ref_positions_less(rhs.ref_positions, lhs.ref_positions)) {
            return false;
        }
        if (lhs.tr_id != rhs.tr_id) {
            return lhs.tr_id < rhs.tr_id;
        }
        if (lhs.depth != rhs.depth) {
            return lhs.depth < rhs.depth;
        }
        return false;
    });

    order.stable_id_by_internal.resize(graph.node_count());
    for (size_t stable_id = 0; stable_id < order.ordered_internal_ids.size(); ++stable_id) {
        order.stable_id_by_internal[order.ordered_internal_ids[stable_id]] = stable_id;
    }

    return order;
}

uint64_t stable_node_id(const StableNodeOrder& order, uint64_t internal_id) {
    return order.stable_id_by_internal[internal_id];
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
        <p class="hint">Expected files: raw.json, clean.json, unitig.json, manifest.json, and optionally flow_paths.json, read_traces.json, locus_traces.json, pair_deletions.json</p>
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
        << "  ],\n"
        << "  \"optional_artifacts\": [\n"
        << "    {\"name\": \"flow_paths\", \"file\": \"flow_paths.json\"},\n"
        << "    {\"name\": \"read_traces\", \"file\": \"read_traces.json\"},\n"
        << "    {\"name\": \"locus_traces\", \"file\": \"locus_traces.json\"},\n"
        << "    {\"name\": \"pair_deletions\", \"file\": \"pair_deletions.json\"}\n"
        << "  ]\n"
        << "}\n";
}

void write_flow_paths_json(const std::string& path,
                           const std::vector<HaplotypePath>& paths) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open flow path JSON: " + path);
    }

    out << "{\n"
        << "  \"graph_kind\": \"flow_paths\",\n"
        << "  \"paths\": [\n";

    for (size_t path_index = 0; path_index < paths.size(); ++path_index) {
        const auto& path_entry = paths[path_index];
        out << "    {\"path_index\": " << path_index
            << ", \"flow\": " << path_entry.flow
            << ", \"unitig_ids\": [";
        for (size_t unitig_index = 0; unitig_index < path_entry.unitig_ids.size(); ++unitig_index) {
            out << path_entry.unitig_ids[unitig_index];
            if (unitig_index + 1 != path_entry.unitig_ids.size()) {
                out << ", ";
            }
        }
        out << "]}";
        out << (path_index + 1 == paths.size() ? "\n" : ",\n");
    }

    out << "  ]\n"
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
                << ", \"ref_positions\": ";
            write_ref_positions_json(out, node.ref_positions);
            out << ", \"tr_id\": " << node.tr_id
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
                << ", \"ref_positions\": ";
            write_ref_positions_json(out, node.ref_positions);
            out << ", \"tr_id\": " << node.tr_id
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

    const StableNodeOrder order = build_stable_node_order(graph);

    out << "H\tVN:Z:1.0\n";

    // S-lines (segments): each node
    for (uint64_t internal_id : order.ordered_internal_ids) {
        const auto& node = graph.node(internal_id);
        out << "S\t" << stable_node_id(order, internal_id) << "\t" << node.kmer
            << "\tDP:f:" << node.depth
            << "\tBB:i:" << (node.is_backbone ? 1 : 0)
            << "\tRP:i:" << node.ref_pos
            << "\tRPS:Z:" << ref_positions_tag(node.ref_positions)
            << "\tTR:i:" << node.tr_id << '\n';
    }

    // L-lines (links): each edge
    std::vector<size_t> edge_order(graph.edge_count());
    std::iota(edge_order.begin(), edge_order.end(), 0);
    std::sort(edge_order.begin(), edge_order.end(), [&graph, &order](size_t lhs_index, size_t rhs_index) {
        const auto& lhs = graph.edges()[lhs_index];
        const auto& rhs = graph.edges()[rhs_index];
        if (stable_node_id(order, lhs.from) != stable_node_id(order, rhs.from)) {
            return stable_node_id(order, lhs.from) < stable_node_id(order, rhs.from);
        }
        if (stable_node_id(order, lhs.to) != stable_node_id(order, rhs.to)) {
            return stable_node_id(order, lhs.to) < stable_node_id(order, rhs.to);
        }
        return lhs.weight < rhs.weight;
    });

    for (size_t edge_index : edge_order) {
        const auto& edge = graph.edges()[edge_index];
        // Overlap = k-1
        int overlap = static_cast<int>(graph.k()) - 1;
        out << "L\t" << stable_node_id(order, edge.from) << "\t+\t"
            << stable_node_id(order, edge.to) << "\t+\t"
            << overlap << "M"
            << "\tRC:i:" << edge.weight << '\n';
    }
}

void write_dbg_json(const std::string& path, const DBG& graph) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open debug JSON: " + path);

    const StableNodeOrder order = build_stable_node_order(graph);

    out << "{\n"
        << "  \"graph_kind\": \"dbg\",\n"
        << "  \"k\": " << graph.k() << ",\n"
        << "  \"nodes\": [\n";

    for (size_t index = 0; index < order.ordered_internal_ids.size(); ++index) {
        uint64_t internal_id = order.ordered_internal_ids[index];
        const auto& node = graph.node(internal_id);
        out << "    {\"id\": " << stable_node_id(order, internal_id)
            << ", \"sequence\": \"" << json_escape(node.kmer) << "\""
            << ", \"ref_pos\": " << node.ref_pos
            << ", \"ref_positions\": ";
        write_ref_positions_json(out, node.ref_positions);
        out << ", \"tr_id\": " << node.tr_id
            << ", \"is_backbone\": " << (node.is_backbone ? "true" : "false")
            << ", \"depth\": " << node.depth
            << ", \"removed\": " << (graph.is_node_removed(internal_id) ? "true" : "false")
            << "}";
        out << (index + 1 == order.ordered_internal_ids.size() ? "\n" : ",\n");
    }

    out << "  ],\n"
        << "  \"edges\": [\n";

    std::vector<size_t> edge_order(graph.edge_count());
    std::iota(edge_order.begin(), edge_order.end(), 0);
    std::sort(edge_order.begin(), edge_order.end(), [&graph, &order](size_t lhs_index, size_t rhs_index) {
        const auto& lhs = graph.edges()[lhs_index];
        const auto& rhs = graph.edges()[rhs_index];
        if (stable_node_id(order, lhs.from) != stable_node_id(order, rhs.from)) {
            return stable_node_id(order, lhs.from) < stable_node_id(order, rhs.from);
        }
        if (stable_node_id(order, lhs.to) != stable_node_id(order, rhs.to)) {
            return stable_node_id(order, lhs.to) < stable_node_id(order, rhs.to);
        }
        return lhs.weight < rhs.weight;
    });
    for (size_t index = 0; index < edge_order.size(); ++index) {
        const auto& edge = graph.edges()[edge_order[index]];
        out << "    {\"from\": " << stable_node_id(order, edge.from)
            << ", \"to\": " << stable_node_id(order, edge.to)
            << ", \"weight\": " << edge.weight
            << "}";
        out << (index + 1 == edge_order.size() ? "\n" : ",\n");
    }

    out << "  ],\n"
        << "  \"haplotype_edges\": [\n";

    std::vector<size_t> hap_edge_order(graph.haplotype_edges().size());
    std::iota(hap_edge_order.begin(), hap_edge_order.end(), 0);
    std::sort(hap_edge_order.begin(), hap_edge_order.end(), [&graph, &order](size_t lhs_index, size_t rhs_index) {
        const auto& lhs = graph.haplotype_edges()[lhs_index];
        const auto& rhs = graph.haplotype_edges()[rhs_index];
        if (stable_node_id(order, lhs.from_node) != stable_node_id(order, rhs.from_node)) {
            return stable_node_id(order, lhs.from_node) < stable_node_id(order, rhs.from_node);
        }
        if (stable_node_id(order, lhs.to_node) != stable_node_id(order, rhs.to_node)) {
            return stable_node_id(order, lhs.to_node) < stable_node_id(order, rhs.to_node);
        }
        return lhs.weight < rhs.weight;
    });
    for (size_t index = 0; index < hap_edge_order.size(); ++index) {
        const auto& edge = graph.haplotype_edges()[hap_edge_order[index]];
        out << "    {\"from\": " << stable_node_id(order, edge.from_node)
            << ", \"to\": " << stable_node_id(order, edge.to_node)
            << ", \"weight\": " << edge.weight
            << "}";
        out << (index + 1 == hap_edge_order.size() ? "\n" : ",\n");
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
            << "\tDP:f:" << u.mean_depth
            << "\tRP:i:" << u.ref_pos
            << "\tRPS:Z:" << ref_positions_tag(u.ref_positions) << '\n';
    }

    for (const auto& e : ug.edges()) {
        out << "L\t" << e.from << "\t+\t" << e.to << "\t+\t0M"
            << "\tRC:i:" << e.weight << '\n';
    }
}

void write_sv_unitig_gfa(const std::string& path, const UnitigGraph& ug) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open GFA: " + path);

    out << "H\tVN:Z:1.0\n";

    for (const auto& u : ug.unitigs()) {
        out << "S\t" << u.id << "\t" << u.sequence
            << "\tDP:f:" << u.mean_depth
            << "\tRP:i:" << u.ref_pos
            << "\tRPS:Z:" << ref_positions_tag(u.ref_positions)
            << "\tSC:Z:" << unitig_support_class_name(u.support_class())
            << "\tBN:i:" << u.backbone_node_count
            << "\tRN:i:" << u.read_node_count
            << "\tCL:z:" << unitig_support_color(u) << '\n';
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
            << ", \"ref_pos\": " << unitig.ref_pos
            << ", \"ref_positions\": ";
        write_ref_positions_json(out, unitig.ref_positions);
        out << ", \"node_ids\": [";
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

void write_sv_unitig_json(const std::string& path, const UnitigGraph& ug) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open debug JSON: " + path);

    out << "{\n"
        << "  \"graph_kind\": \"unitig_sv\",\n"
        << "  \"unitigs\": [\n";

    const auto& unitigs = ug.unitigs();
    for (size_t index = 0; index < unitigs.size(); ++index) {
        const auto& unitig = unitigs[index];
        out << "    {\"id\": " << unitig.id
            << ", \"sequence\": \"" << json_escape(unitig.sequence) << "\""
            << ", \"mean_depth\": " << unitig.mean_depth
            << ", \"ref_pos\": " << unitig.ref_pos
            << ", \"ref_positions\": ";
        write_ref_positions_json(out, unitig.ref_positions);
        out << ", \"node_ids\": [";
        for (size_t node_index = 0; node_index < unitig.node_ids.size(); ++node_index) {
            out << unitig.node_ids[node_index];
            if (node_index + 1 != unitig.node_ids.size()) {
                out << ", ";
            }
        }
        out << "]"
            << ", \"support_class\": \"" << unitig_support_class_name(unitig.support_class()) << "\""
            << ", \"backbone_node_count\": " << unitig.backbone_node_count
            << ", \"read_node_count\": " << unitig.read_node_count
            << ", \"color\": \"" << unitig_support_color(unitig) << "\""
            << "}";
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
            if (node.has_ref_pos_in_range(start, end)) {
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
                "",
                {node.ref_positions.begin(), node.ref_positions.end()}
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

void write_flow_path_artifacts(const DebugArtifactsConfig& config,
                               const std::vector<HaplotypePath>& paths) {
    if (!config.should_write() || paths.empty()) {
        return;
    }

    ensure_debug_output_dir(config);
    write_flow_paths_json(flow_paths_json_path(config), paths);
}

void write_vcf(const std::string& path,
               const std::vector<StructuralVariantCall>& calls,
               const std::string& source) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open VCF: " + path);
    }

    out << "##fileformat=VCFv4.3\n";
    out << "##source=" << escape_vcf_field(source) << "\n";
    out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"1-based inclusive end position of the reference allele\">\n";
    out << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"ALT length minus REF length\">\n";
    out << "##INFO=<ID=SUPPORT,Number=1,Type=Float,Description=\"Representative alternate-path support score (minimum mean depth across non-backbone unitigs)\">\n";
    out << "##INFO=<ID=SRC_UID,Number=1,Type=Integer,Description=\"Source backbone unitig ID for the reference interval\">\n";
    out << "##INFO=<ID=SNK_UID,Number=1,Type=Integer,Description=\"Sink backbone unitig ID for the reference interval\">\n";
    out << "##INFO=<ID=SRC_REF_POS,Number=1,Type=Integer,Description=\"1-based reference anchor coordinate of the source backbone unitig\">\n";
    out << "##INFO=<ID=SNK_REF_POS,Number=1,Type=Integer,Description=\"1-based reference anchor coordinate of the sink backbone unitig\">\n";
    out << "##INFO=<ID=CALL_SOURCE,Number=1,Type=String,Description=\"Primary evidence source for the call (PAIR, PATH_PAIR)\">\n";
    out << "##INFO=<ID=PAIR_SUPPORT,Number=1,Type=Integer,Description=\"Number of abnormal read pairs supporting a pair-rescued deletion\">\n";
    out << "##INFO=<ID=PAIR_MAX_TLEN,Number=1,Type=Integer,Description=\"Maximum observed template length among supporting read pairs\">\n";
    out << "##INFO=<ID=PAIR_ANCHOR_GAP,Number=1,Type=Integer,Description=\"Maximum observed mapped reference gap between paired-read anchors\">\n";
    out << "##INFO=<ID=PAIR_IMPLIED_DEL,Number=1,Type=Integer,Description=\"Maximum deletion span implied by abnormal template length relative to the local fragment model\">\n";
    out << "##INFO=<ID=PAIR_COV_RATIO,Number=1,Type=Float,Description=\"Interior-to-flank depth ratio for the pair-supported reference interval\">\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    std::vector<StructuralVariantCall> sorted_calls = calls;
    std::sort(sorted_calls.begin(), sorted_calls.end(),
              [](const StructuralVariantCall& lhs, const StructuralVariantCall& rhs) {
                  if (lhs.chrom != rhs.chrom) {
                      return lhs.chrom < rhs.chrom;
                  }
                  if (lhs.pos != rhs.pos) {
                      return lhs.pos < rhs.pos;
                  }
                  return lhs.end < rhs.end;
              });

    for (const auto& call : sorted_calls) {
        out << escape_vcf_field(call.chrom) << '\t'
            << std::max<int32_t>(1, call.pos + 1) << '\t'
            << escape_vcf_field(call.id) << '\t'
            << escape_vcf_field(call.ref) << '\t'
            << escape_vcf_field(call.alt) << '\t'
            << ".\t"
            << escape_vcf_field(call.filter) << '\t'
            << build_vcf_info_field(call) << '\n';
    }
}

} // namespace sharda
