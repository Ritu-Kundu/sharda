#pragma once

#include "graph/types.h"
#include <string>
#include <vector>

namespace sharda {

class UnitigGraph;

std::vector<StructuralVariantCall> call_structural_variants(
    const UnitigGraph& graph,
    const std::string& chrom,
    int32_t coord_offset = 0,
    size_t max_paths_per_interval = 256);

std::vector<StructuralVariantCall> call_pair_supported_deletions(
    const UnitigGraph& graph,
    const std::vector<ReadPair>& read_pairs,
    const std::string& chrom,
    int32_t coord_offset = 0,
    size_t min_supporting_pairs = 2,
    double fragment_mean_fallback = -1.0,
    double fragment_sd_fallback = -1.0);

std::vector<StructuralVariantCall> merge_structural_variant_sources(
    std::vector<StructuralVariantCall> sequence_calls,
    const std::vector<StructuralVariantCall>& pair_calls);

std::vector<StructuralVariantCall> collapse_structural_variant_calls(
    const std::vector<StructuralVariantCall>& calls);

} // namespace sharda