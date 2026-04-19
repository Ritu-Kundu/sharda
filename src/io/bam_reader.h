#pragma once

#include "graph/types.h"
#include <functional>
#include <string>
#include <cstdint>

namespace sharda {

/// Iterate over read pairs in a name-sorted BAM.
/// Calls `callback` for each read pair.  Skips secondary/supplementary records.
/// Throws on I/O errors.
void iterate_read_pairs(const std::string& bam_path,
                        std::function<void(ReadPair&&)> callback);

/// Extract reads overlapping a genomic region from an indexed BAM,
/// name-sort them, and write to a temporary name-sorted BAM.
/// Requires a coordinate-sorted BAM with an associated .bai index.
/// region string format: "chr1:1000-2000" (1-based, inclusive).
void create_region_bam(const std::string& src_bam,
                       const std::string& chrom,
                       int32_t start, int32_t end,
                       const std::string& dest_bam);

} // namespace sharda
