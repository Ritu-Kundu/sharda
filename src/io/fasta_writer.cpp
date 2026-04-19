#include "io/fasta_writer.h"
#include <fstream>
#include <stdexcept>

namespace sharda {

void write_fasta(const std::string& path,
                 const std::vector<std::pair<std::string, std::string>>& entries) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open output FASTA: " + path);

    for (const auto& [name, seq] : entries) {
        out << '>' << name << '\n';
        // Wrap at 80 characters
        for (size_t i = 0; i < seq.size(); i += 80) {
            out << seq.substr(i, 80) << '\n';
        }
    }
}

} // namespace sharda
