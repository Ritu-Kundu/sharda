#pragma once

#include "graph/types.h"
#include <string>
#include <vector>
#include <utility>

namespace sharda {

/// Parameters for assembling a single region.
struct RegionParams {
    std::string region_name;  // e.g. "chr1:10000-20000"
    std::string ref_seq;      // extracted reference subsequence (local coords)
    std::string bam_path;     // name-sorted BAM for this region
    int32_t     coord_offset; // genomic start of the extracted region
    std::vector<TandemRepeat> trs; // TRs in local coordinates
    int         ploidy  = 2;
    int         k       = 121;
    bool        debug   = false;
    std::string debug_dir;    // if non-empty, write GFA files here
};

/// Result from assembling a single region.
struct RegionResult {
    std::string region_name;
    std::vector<std::pair<std::string, std::string>> haplotypes; // (name, sequence)
    bool        success = false;
    std::string error;
};

/// Assemble haplotypes for a single target region.
/// Thread-safe: each invocation operates on independent data.
RegionResult assemble_region(const RegionParams& params);

} // namespace sharda
