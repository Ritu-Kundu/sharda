#include "assembly/read_adder.h"
#include "assembly/read_classifier.h"
#include "assembly/anchor_chain.h"
#include "util/kmer.h"
#include <spdlog/spdlog.h>
#include <algorithm>

namespace sharda {

namespace {

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
    trace_record->raw_nodes.push_back({
        node.id,
        node.kmer,
        created,
        node.is_backbone,
        node.ref_pos,
        node.tr_id,
        false,
        UINT64_MAX,
        ""
    });
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

/// Walk the CIGAR to produce (read_offset, ref_pos) pairs for each aligned base.
/// Includes soft-clipped bases with ref_pos = -1.
struct AlignedBase {
    int read_offset;
    int32_t ref_pos; // -1 for soft-clip/insertion
};

std::vector<AlignedBase> cigar_walk(const AlignedRead& read) {
    std::vector<AlignedBase> bases;
    bases.reserve(read.seq.size());
    int roff = 0;       // read offset
    int32_t rpos = read.ref_start; // ref position

    for (const auto& c : read.cigar) {
        switch (c.op) {
        case CigarOp::M:
        case CigarOp::EQ:
        case CigarOp::X:
            for (uint32_t i = 0; i < c.len; ++i) {
                bases.push_back({roff++, rpos++});
            }
            break;
        case CigarOp::I:
            for (uint32_t i = 0; i < c.len; ++i) {
                bases.push_back({roff++, -1});
            }
            break;
        case CigarOp::D:
        case CigarOp::N:
            rpos += c.len;
            break;
        case CigarOp::S:
            for (uint32_t i = 0; i < c.len; ++i) {
                bases.push_back({roff++, -1});
            }
            break;
        case CigarOp::H:
        case CigarOp::P:
            break; // hard-clip / padding: skip
        }
    }
    return bases;
}

/// Add an ORR-style read path using alignment positions.
/// Returns (first_node_id, last_node_id) added to the path.
std::pair<uint64_t, uint64_t> add_orr_path(
    const AlignedRead& read,
    DBG& graph,
    const std::vector<AlignedBase>& aligned_bases,
    ReadTraceRecord* trace_record = nullptr)
{
    int k = graph.k();
    auto kmers = extract_kmers(read.seq, k);
    if (kmers.empty()) return {UINT64_MAX, UINT64_MAX};

    // For each kmer, determine its ref_pos from the first base of the kmer
    // that has a valid ref_pos (the aligned start of the kmer).
    uint64_t first_nid = UINT64_MAX;
    uint64_t prev_nid  = UINT64_MAX;

    for (int ki = 0; ki < static_cast<int>(kmers.size()); ++ki) {
        // Find the ref_pos for this kmer: use the ref_pos of aligned_bases[ki]
        // if available
        int32_t kmer_ref_pos = -1;
        if (ki < static_cast<int>(aligned_bases.size())) {
            kmer_ref_pos = aligned_bases[ki].ref_pos;
        }

        uint64_t nid;
        bool created = false;
        if (kmer_ref_pos >= 0) {
            // Try to match to a backbone node
            nid = graph.backbone_node_at(kmer_ref_pos);
            if (nid != UINT64_MAX && graph.node(nid).kmer == kmers[ki]) {
                // Reuse backbone node
                graph.node_mut(nid).depth++;
            } else {
                // kmer differs from backbone at this position: novel node
                created = graph.find_read_node(kmers[ki]) == UINT64_MAX;
                nid = graph.add_read_node(kmers[ki]);
                graph.node_mut(nid).depth++;
            }
        } else {
            // Soft-clip or insertion: novel node
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

/// Add an IRR-style read path using anchor chaining.
std::pair<uint64_t, uint64_t> add_irr_path(
    const AlignedRead& read,
    int tr_id,
    DBG& graph,
    const std::vector<AlignedBase>& aligned_bases,
    ReadTraceRecord* trace_record = nullptr)
{
    int k = graph.k();
    auto kmers = extract_kmers(read.seq, k);
    if (kmers.empty()) return {UINT64_MAX, UINT64_MAX};

    auto anchors = find_and_chain_anchors(kmers, graph, tr_id);

    if (anchors.empty()) {
        // Fall back to ORR-style
        return add_orr_path(read, graph, aligned_bases, trace_record);
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

    auto ab1 = cigar_walk(pair.read1);
    auto ab2 = cigar_walk(pair.read2);

    std::pair<uint64_t, uint64_t> path1, path2;
    ReadTraceRecord* trace1 = start_trace_record(trace_sink, pair.read1, "read1", cls1);
    ReadTraceRecord* trace2 = start_trace_record(trace_sink, pair.read2, "read2", cls2);

    if (cls1.type == ReadType::IRR) {
        path1 = add_irr_path(pair.read1, cls1.tr_id, graph, ab1, trace1);
    } else {
        path1 = add_orr_path(pair.read1, graph, ab1, trace1);
    }

    if (cls2.type == ReadType::IRR) {
        path2 = add_irr_path(pair.read2, cls2.tr_id, graph, ab2, trace2);
    } else {
        path2 = add_orr_path(pair.read2, graph, ab2, trace2);
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
