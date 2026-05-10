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

std::vector<StructuralVariantCall> collapse_structural_variant_calls(
    const std::vector<StructuralVariantCall>& calls);

} // namespace sharda