#pragma once

#include "graph/dbg.h"

namespace sharda {

/// Clean the graph by removing tips, low-weight edges, and error bubbles.
/// Operates iteratively until no more changes.
/// mean_read_length: used for tip removal threshold.
void clean_graph(DBG& graph, int mean_read_length);

} // namespace sharda
