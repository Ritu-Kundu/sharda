#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace sharda {

// ── Genomic region ──────────────────────────────────────────────────────────
struct TargetRegion {
    std::string chrom;
    int32_t start = 0; // 0-based, inclusive
    int32_t end   = 0; // 0-based, exclusive
};

// ── Tandem repeat from BED ──────────────────────────────────────────────────
struct TandemRepeat {
    std::string chrom;
    int32_t start = 0; // 0-based, inclusive
    int32_t end   = 0; // 0-based, exclusive
    int         id    = -1;
};

// ── Graph node ──────────────────────────────────────────────────────────────
struct Node {
    uint64_t    id        = 0;
    std::string kmer;
    int32_t     ref_pos   = -1; // -1 for non-backbone
    int         tr_id     = -1; // -1 if not in any TR
    bool        is_backbone = false;
    uint32_t    depth     = 0;  // number of reads covering this node
};

// ── Directed edge in the DBG ────────────────────────────────────────────────
struct Edge {
    uint64_t from   = 0;
    uint64_t to     = 0;
    uint32_t weight = 0;
};

// ── Haplotype edge (read-pair phasing) ──────────────────────────────────────
struct HaplotypeEdge {
    uint64_t from_node = 0;
    uint64_t to_node   = 0;
    uint32_t weight    = 0;
};

// ── Aligned read ────────────────────────────────────────────────────────────
enum class CigarOp : uint8_t {
    M = 0, I = 1, D = 2, N = 3, S = 4, H = 5, P = 6, EQ = 7, X = 8
};

struct CigarElement {
    CigarOp  op;
    uint32_t len;
};

struct AlignedRead {
    std::string              name;
    std::string              seq;
    std::string              qual;
    std::vector<CigarElement> cigar;
    int32_t                  ref_start  = 0; // 0-based leftmost mapping pos
    int32_t                  ref_end    = 0; // 0-based exclusive
    uint16_t                 flag       = 0;
    bool                     has_sa_tag = false; // supplementary alignment tag

    bool is_unmapped()      const { return flag & 0x4; }
    bool is_reverse()       const { return flag & 0x10; }
    bool is_secondary()     const { return flag & 0x100; }
    bool is_supplementary() const { return flag & 0x800; }
    bool is_proper_pair()   const { return flag & 0x2; }
    bool is_read1()         const { return flag & 0x40; }
};

// ── Read pair ───────────────────────────────────────────────────────────────
struct ReadPair {
    AlignedRead read1;
    AlignedRead read2;
};

// ── Read classification ─────────────────────────────────────────────────────
enum class ReadType { ORR, IRR };

struct ReadClassification {
    ReadType type   = ReadType::ORR;
    int      tr_id  = -1; // TR id for IRR reads
    bool     is_evidence = false;
};

// ── Unitig ──────────────────────────────────────────────────────────────────
struct Unitig {
    uint64_t              id = 0;
    std::string           sequence;
    double                mean_depth = 0.0;
    std::vector<uint64_t> node_ids; // constituent node IDs
};

} // namespace sharda
