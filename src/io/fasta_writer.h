#pragma once

#include <string>
#include <vector>
#include <utility>

namespace sharda {

/// Write sequences to a FASTA file. Each entry is (name, sequence).
void write_fasta(const std::string& path,
                 const std::vector<std::pair<std::string, std::string>>& entries);

} // namespace sharda
