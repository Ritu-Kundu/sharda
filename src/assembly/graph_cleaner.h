#pragma once

#include "graph/dbg.h"

namespace sharda {

struct GraphCleaningOptions {
	bool preserve_backbone_edges = false;
};

/// Clean the graph by removing tips, low-weight edges, and error bubbles.
/// Operates iteratively until no more changes.
/// mean_read_length: used for tip removal threshold.
void clean_graph(DBG& graph,
				 int mean_read_length,
				 GraphCleaningOptions options = {});

} // namespace sharda
