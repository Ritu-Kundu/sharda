#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace sharda {

/// Extract all k-mers (forward only) from a sequence.
/// Returns empty vector if seq.size() < k.
std::vector<std::string> extract_kmers(std::string_view seq, int k);

} // namespace sharda
