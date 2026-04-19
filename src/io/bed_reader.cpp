#include "io/bed_reader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace sharda {

std::vector<TandemRepeat> read_bed(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open BED: " + path);

    std::vector<TandemRepeat> trs;
    std::string line;
    int next_id = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        TandemRepeat tr;
        if (!(iss >> tr.chrom >> tr.start >> tr.end))
            throw std::runtime_error("Malformed BED line: " + line);
        tr.id = next_id++;
        trs.push_back(std::move(tr));
    }
    return trs;
}

std::vector<TargetRegion> read_target_regions(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open target BED: " + path);

    std::vector<TargetRegion> regions;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        TargetRegion r;
        if (!(iss >> r.chrom >> r.start >> r.end))
            throw std::runtime_error("Malformed BED line: " + line);
        regions.push_back(std::move(r));
    }
    return regions;
}

std::vector<TandemRepeat> filter_trs_for_region(
    const std::vector<TandemRepeat>& all_trs,
    const TargetRegion& region,
    int32_t padding)
{
    int32_t ext_start = std::max(static_cast<int32_t>(0), region.start - padding);
    int32_t ext_end   = region.end + padding;

    std::vector<TandemRepeat> local;
    int next_id = 0;
    for (const auto& tr : all_trs) {
        if (tr.chrom != region.chrom) continue;
        if (tr.start >= ext_end || tr.end <= ext_start) continue; // no overlap

        TandemRepeat adj = tr;
        adj.start = std::max(tr.start, ext_start) - ext_start;
        adj.end   = std::min(tr.end,   ext_end)   - ext_start;
        adj.id    = next_id++;
        local.push_back(std::move(adj));
    }
    return local;
}

} // namespace sharda
