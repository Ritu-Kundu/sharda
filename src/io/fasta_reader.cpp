#include "io/fasta_reader.h"
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <htslib/faidx.h>

namespace sharda {

std::pair<std::string, std::string> read_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open FASTA: " + path);

    std::string name, line, seq;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!name.empty()) break; // only first record
            // name is everything after '>' up to first whitespace
            size_t end = line.find_first_of(" \t", 1);
            name = line.substr(1, end - 1);
        } else {
            // Remove whitespace and append
            line.erase(std::remove_if(line.begin(), line.end(),
                        [](unsigned char c){ return std::isspace(c); }),
                        line.end());
            seq += line;
        }
    }
    if (name.empty() || seq.empty())
        throw std::runtime_error("No sequence found in FASTA: " + path);

    // Upper-case normalise
    std::transform(seq.begin(), seq.end(), seq.begin(),
                   [](unsigned char c){ return std::toupper(c); });
    return {name, seq};
}

std::string read_fasta_region(const std::string& fasta_path,
                              const std::string& chrom,
                              int32_t start, int32_t end) {
    faidx_t* fai = fai_load(fasta_path.c_str());
    if (!fai)
        throw std::runtime_error("Cannot load FASTA index for: " + fasta_path
                                 + " (is the .fai file present?)");

    int len = 0;
    // faidx_fetch_seq uses 0-based inclusive start, inclusive end
    char* seq = faidx_fetch_seq(fai, chrom.c_str(), start, end - 1, &len);
    if (!seq || len <= 0) {
        fai_destroy(fai);
        throw std::runtime_error("Cannot fetch region " + chrom + ":"
                                 + std::to_string(start) + "-" + std::to_string(end)
                                 + " from " + fasta_path);
    }

    std::string result(seq, static_cast<size_t>(len));
    free(seq);
    fai_destroy(fai);

    // Upper-case normalise
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c){ return std::toupper(c); });
    return result;
}

} // namespace sharda
