#include "assembly/read_adder.h"
#include "assembly/read_classifier.h"
#include "assembly/anchor_chain.h"
#include "util/kmer.h"
#include <spdlog/spdlog.h>
#include <algorithm>

namespace sharda {

namespace {

std::vector<int32_t> to_ref_position_vector(const std::set<int32_t>& ref_positions) {
    return {ref_positions.begin(), ref_positions.end()};
}

bool should_trace_read(const ReadTraceSink& trace_sink, const std::string& read_name) {
    if (!trace_sink.enabled()) {
        return false;
    }
    return std::find(trace_sink.traced_reads->begin(),
                     trace_sink.traced_reads->end(),
                     read_name) != trace_sink.traced_reads->end();
}

std::string read_type_name(ReadType type) {
    return type == ReadType::IRR ? "IRR" : "ORR";
}

void append_trace_node(ReadTraceRecord* trace_record,
                       const DBG& graph,
                       uint64_t node_id,
                       bool created) {
    if (trace_record == nullptr || node_id == UINT64_MAX) {
        return;
    }

    const auto& node = graph.node(node_id);

    ReadTraceNode trace_node;
    trace_node.node_id = node.id;
    trace_node.sequence = node.kmer;
    trace_node.created = created;
    trace_node.is_backbone = node.is_backbone;
    trace_node.ref_pos = node.ref_pos;
    trace_node.tr_id = node.tr_id;
    trace_node.removed_after_clean = false;
    trace_node.unitig_id = UINT64_MAX;
    trace_node.unitig_sequence.clear();
    trace_node.ref_positions = to_ref_position_vector(node.ref_positions);
    trace_record->raw_nodes.push_back(std::move(trace_node));
}

ReadTraceRecord* start_trace_record(const ReadTraceSink& trace_sink,
                                    const AlignedRead& read,
                                    const char* mate_label,
                                    const ReadClassification& classification) {
    if (!should_trace_read(trace_sink, read.name)) {
        return nullptr;
    }

    trace_sink.records->push_back({
        read.name,
        mate_label,
        read_type_name(classification.type),
        classification.is_evidence,
        classification.tr_id,
        {}
    });
    return &trace_sink.records->back();
}

uint64_t choose_path_node_for_kmer(const std::string& kmer,
                                   int32_t implied_ref_pos,
                                   DBG& graph,
                                   bool& created) {
    uint64_t backbone_nid = graph.closest_backbone_node_for_kmer(kmer, implied_ref_pos);
    if (backbone_nid != UINT64_MAX) {
        created = false;
        return backbone_nid;
    }

    uint64_t existing_read_nid = graph.find_read_node(kmer);
    created = existing_read_nid == UINT64_MAX;
    uint64_t nid = graph.add_read_node(kmer);
    graph.add_node_ref_pos(nid, implied_ref_pos);
    return nid;
}

/// Add an ORR-style read path using implied coordinates from the alignment start.
/// Returns (first_node_id, last_node_id) added to the path.
std::pair<uint64_t, uint64_t> add_orr_path(
    const AlignedRead& read,
    DBG& graph,
    ReadTraceRecord* trace_record = nullptr)
{
    int k = graph.k();
    auto kmers = extract_kmers(read.seq, k);
    if (kmers.empty()) return {UINT64_MAX, UINT64_MAX};

    uint64_t first_nid = UINT64_MAX;
    uint64_t prev_nid  = UINT64_MAX;

    for (int ki = 0; ki < static_cast<int>(kmers.size()); ++ki) {
        int32_t implied_ref_pos = read.ref_start + ki;
        bool created = false;
        uint64_t nid = choose_path_node_for_kmer(kmers[ki], implied_ref_pos, graph, created);
        graph.node_mut(nid).depth++;

        if (first_nid == UINT64_MAX) {
            first_nid = nid;

            if (!graph.node(nid).is_backbone && implied_ref_pos > 0) {
                uint64_t predecessor = graph.backbone_node_at(implied_ref_pos - 1);
                if (predecessor != UINT64_MAX && predecessor != nid) {
                    graph.add_edge(predecessor, nid);
                }
            }
        }

        append_trace_node(trace_record, graph, nid, created);

        if (prev_nid != UINT64_MAX && prev_nid != nid) {
            graph.add_edge(prev_nid, nid);
        }
        prev_nid = nid;
    }

    return {first_nid, prev_nid};
}

/// Add an IRR-style read path using anchor chaining.
std::pair<uint64_t, uint64_t> add_irr_path(
    const AlignedRead& read,
    int tr_id,
    DBG& graph,
    ReadTraceRecord* trace_record = nullptr)
{
    int k = graph.k();
    auto kmers = extract_kmers(read.seq, k);
    if (kmers.empty()) return {UINT64_MAX, UINT64_MAX};

    auto anchors = find_and_chain_anchors(kmers, graph, tr_id);

    if (anchors.empty()) {
        // Fall back to ORR-style implied-coordinate placement.
        return add_orr_path(read, graph, trace_record);
    }

    uint64_t first_nid = UINT64_MAX;
    uint64_t prev_nid  = UINT64_MAX;

    // Walk through read kmers, using anchors where available
    // Between anchors, add in traditional DBG style (hash-based)
    int ai = 0; // anchor index
    for (int ki = 0; ki < static_cast<int>(kmers.size()); ++ki) {
        uint64_t nid;
        bool created = false;

        if (ai < static_cast<int>(anchors.size()) && anchors[ai].first == ki) {
            // This kmer is an anchor — use the backbone node
            nid = anchors[ai].second;
            graph.node_mut(nid).depth++;
            ++ai;
        } else {
            // Between anchors: add as read node (DBG-style)
            created = graph.find_read_node(kmers[ki]) == UINT64_MAX;
            nid = graph.add_read_node(kmers[ki]);
            graph.node_mut(nid).depth++;
        }

        append_trace_node(trace_record, graph, nid, created);

        if (first_nid == UINT64_MAX) first_nid = nid;
        if (prev_nid != UINT64_MAX && prev_nid != nid) {
            graph.add_edge(prev_nid, nid);
        }
        prev_nid = nid;
    }

    return {first_nid, prev_nid};
}

} // anonymous namespace

void add_read_pair(ReadPair& pair, DBG& graph,
                   const std::vector<TandemRepeat>& trs) {
    add_read_pair(pair, graph, trs, {});
}

void add_read_pair(ReadPair& pair, DBG& graph,
                   const std::vector<TandemRepeat>& trs,
                   const ReadTraceSink& trace_sink) {
    auto cls1 = classify_read(pair.read1, trs);
    auto cls2 = classify_read(pair.read2, trs);

    if (trace_sink.enabled()) {
        size_t extra_records = 0;
        extra_records += should_trace_read(trace_sink, pair.read1.name) ? 1u : 0u;
        extra_records += should_trace_read(trace_sink, pair.read2.name) ? 1u : 0u;
        if (extra_records > 0) {
            trace_sink.records->reserve(trace_sink.records->size() + extra_records);
        }
    }

    std::pair<uint64_t, uint64_t> path1, path2;
    ReadTraceRecord* trace1 = start_trace_record(trace_sink, pair.read1, "read1", cls1);
    ReadTraceRecord* trace2 = start_trace_record(trace_sink, pair.read2, "read2", cls2);

    if (cls1.type == ReadType::IRR) {
        path1 = add_irr_path(pair.read1, cls1.tr_id, graph, trace1);
    } else {
        path1 = add_orr_path(pair.read1, graph, trace1);
    }

    if (cls2.type == ReadType::IRR) {
        path2 = add_irr_path(pair.read2, cls2.tr_id, graph, trace2);
    } else {
        path2 = add_orr_path(pair.read2, graph, trace2);
    }

    // Add haplotype edge if at least one read is evidence
    if (cls1.is_evidence || cls2.is_evidence) {
        // Determine left/right by alignment position
        uint64_t left_last, right_first;
        if (pair.read1.ref_start <= pair.read2.ref_start) {
            left_last   = path1.second; // last node of left read
            right_first = path2.first;  // first node of right read
        } else {
            left_last   = path2.second;
            right_first = path1.first;
        }
        if (left_last != UINT64_MAX && right_first != UINT64_MAX) {
            graph.add_haplotype_edge(left_last, right_first);
        }
    }
}

} // namespace sharda
