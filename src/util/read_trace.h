#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace sharda {

struct ReadTraceNode {
    uint64_t    node_id = UINT64_MAX;
    std::string sequence;
    bool        created = false;
    bool        is_backbone = false;
    int32_t     ref_pos = -1;
    int         tr_id = -1;
    bool        removed_after_clean = false;
    uint64_t    unitig_id = UINT64_MAX;
    std::string unitig_sequence;
    std::vector<int32_t> ref_positions;
};

struct ReadTraceRecord {
    std::string read_name;
    std::string mate_label;
    std::string read_type;
    bool        is_evidence = false;
    int         tr_id = -1;
    std::vector<ReadTraceNode> raw_nodes;
};

struct ReadTraceSink {
    const std::vector<std::string>* traced_reads = nullptr;
    std::vector<ReadTraceRecord>* records = nullptr;

    bool enabled() const {
        return traced_reads != nullptr && records != nullptr;
    }
};

} // namespace sharda