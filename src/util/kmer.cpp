#include "util/kmer.h"

namespace sharda {

std::vector<std::string> extract_kmers(std::string_view seq, int k) {
    std::vector<std::string> kmers;
    if (k <= 0 || static_cast<int>(seq.size()) < k) return kmers;
    kmers.reserve(seq.size() - k + 1);
    for (size_t i = 0; i + k <= seq.size(); ++i) {
        kmers.emplace_back(seq.substr(i, k));
    }
    return kmers;
}

} // namespace sharda
