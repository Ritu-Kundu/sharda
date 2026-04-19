#pragma once

#include "graph/types.h"
#include <string>
#include <vector>

namespace sharda {

/// Read tandem repeat annotations from a BED file.
/// Assigns sequential IDs starting from 0.
std::vector<TandemRepeat> read_bed(const std::string& path);

/// Read target regions from a BED file.
std::vector<TargetRegion> read_target_regions(const std::string& path);

/// Filter TRs that overlap a padded target region and adjust to local coords.
/// Returns TRs with coordinates relative to (region.start - padding).
/// The TR ids are reassigned sequentially starting from 0.
std::vector<TandemRepeat> filter_trs_for_region(
    const std::vector<TandemRepeat>& all_trs,
    const TargetRegion& region,
    int32_t padding);

} // namespace sharda
