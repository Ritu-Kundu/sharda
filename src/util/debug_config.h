#pragma once

#include "util/locus_trace.h"
#include <string>
#include <vector>

namespace sharda {

struct DebugArtifactsConfig {
    bool        enabled        = false;
    std::string output_dir;
    bool        emit_gfa       = true;
    bool        emit_json      = true;
    bool        emit_html_view = true;
    std::vector<std::string> traced_reads;
    std::vector<LocusTraceRequest> traced_loci;

    bool should_write() const {
        return enabled && !output_dir.empty();
    }

    bool should_trace_reads() const {
        return should_write() && !traced_reads.empty();
    }

    bool should_trace_loci() const {
        return should_write() && !traced_loci.empty();
    }
};

} // namespace sharda