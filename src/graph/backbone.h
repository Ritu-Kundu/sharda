#pragma once

#include "graph/dbg.h"
#include "graph/types.h"
#include <string>
#include <vector>

namespace sharda {

/// Build the backbone (linear reference path) in the DBG.
/// Creates one node per ref position and connects them sequentially.
void build_backbone(DBG& graph,
                    const std::string& ref_seq,
                    const std::vector<TandemRepeat>& trs);

} // namespace sharda
