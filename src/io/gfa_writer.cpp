#include "io/gfa_writer.h"
#include "graph/dbg.h"
#include "graph/unitig_graph.h"
#include <fstream>
#include <stdexcept>

namespace sharda {

void write_gfa(const std::string& path, const DBG& graph) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open GFA: " + path);

    out << "H\tVN:Z:1.0\n";

    // S-lines (segments): each node
    for (const auto& node : graph.nodes()) {
        out << "S\t" << node.id << "\t" << node.kmer
            << "\tDP:f:" << node.depth
            << "\tBB:i:" << (node.is_backbone ? 1 : 0)
            << "\tRP:i:" << node.ref_pos
            << "\tTR:i:" << node.tr_id << '\n';
    }

    // L-lines (links): each edge
    for (const auto& edge : graph.edges()) {
        // Overlap = k-1
        int overlap = static_cast<int>(graph.k()) - 1;
        out << "L\t" << edge.from << "\t+\t" << edge.to << "\t+\t"
            << overlap << "M"
            << "\tRC:i:" << edge.weight << '\n';
    }
}

void write_unitig_gfa(const std::string& path, const UnitigGraph& ug) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open GFA: " + path);

    out << "H\tVN:Z:1.0\n";

    for (const auto& u : ug.unitigs()) {
        out << "S\t" << u.id << "\t" << u.sequence
            << "\tDP:f:" << u.mean_depth << '\n';
    }

    for (const auto& e : ug.edges()) {
        out << "L\t" << e.from << "\t+\t" << e.to << "\t+\t0M"
            << "\tRC:i:" << e.weight << '\n';
    }
}

} // namespace sharda
