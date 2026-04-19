#pragma once

#include <cstdint>
#include <string>
#include <utility>

namespace sharda {

/// Read the first sequence from a FASTA file.
/// Returns (name, sequence). Throws on error.
std::pair<std::string, std::string> read_fasta(const std::string& path);

/// Extract a region from an indexed FASTA (requires .fai index).
/// chrom: chromosome name, start/end: 0-based half-open [start, end).
/// Returns the upper-cased subsequence. Throws on error.
std::string read_fasta_region(const std::string& fasta_path,
                              const std::string& chrom,
                              int32_t start, int32_t end);

} // namespace sharda
