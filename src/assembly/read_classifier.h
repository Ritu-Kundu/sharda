#pragma once

#include "graph/types.h"
#include <vector>

namespace sharda {

/// Classify a single read as evidence or not, and as ORR or IRR.
ReadClassification classify_read(const AlignedRead& read,
                                 const std::vector<TandemRepeat>& trs);

} // namespace sharda
