#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace sharda {

struct LocusTraceRequest {
    int32_t start = 0;
    int32_t length = 0;
};

struct LocusTraceNode {
    uint64_t    node_id = UINT64_MAX;
    std::string sequence;
    bool        is_backbone = false;
    int32_t     ref_pos = -1;
    int         tr_id = -1;
    bool        removed_after_clean = false;
    uint64_t    unitig_id = UINT64_MAX;
    std::string unitig_sequence;
    std::vector<int32_t> ref_positions;
};

struct LocusTraceRecord {
    int32_t                    local_start = 0;
    int32_t                    length = 0;
    int32_t                    global_start = 0;
    std::string                reference_sequence;
    std::vector<LocusTraceNode> raw_nodes;
};

} // namespace sharda