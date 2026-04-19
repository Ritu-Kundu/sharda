#pragma once

#include "graph/dbg.h"
#include "graph/types.h"
#include <vector>

namespace sharda {

/// Add a read pair to the graph. Classifies reads, adds nodes/edges,
/// and adds a haplotype edge if appropriate.
void add_read_pair(ReadPair& pair, DBG& graph,
                   const std::vector<TandemRepeat>& trs);

} // namespace sharda
