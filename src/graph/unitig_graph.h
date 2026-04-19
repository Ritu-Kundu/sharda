#pragma once

#include "graph/types.h"
#include "graph/dbg.h"
#include <vector>
#include <unordered_map>

namespace sharda {

/// Compacted unitig graph built from a cleaned DBG.
class UnitigGraph {
public:
    /// Build from a cleaned DBG. Returns false if cycles are detected.
    bool build(const DBG& source);

    const std::vector<Unitig>& unitigs() const { return unitigs_; }
    const std::vector<Edge>& edges() const { return edges_; }
    const std::vector<HaplotypeEdge>& haplotype_edges() const { return hap_edges_; }

    /// Check if the graph has cycles (set during build).
    bool has_cycles() const { return has_cycles_; }

    /// Map from original DBG node id to unitig id.
    uint64_t node_to_unitig(uint64_t node_id) const;

    size_t unitig_count() const { return unitigs_.size(); }

    const Unitig& unitig(uint64_t id) const { return unitigs_[id]; }

private:
    std::vector<Unitig> unitigs_;
    std::vector<Edge>   edges_;
    std::vector<HaplotypeEdge> hap_edges_;
    std::unordered_map<uint64_t, uint64_t> node_to_unitig_;
    bool has_cycles_ = false;

    bool detect_cycles() const;
};

} // namespace sharda
