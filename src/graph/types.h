#pragma once

#include <cstdint>
#include <set>
#include <string>
#include <vector>

namespace sharda {

enum class ExecutionMode : uint8_t {
    Haplotype,
    Sv,
    Both,
};

enum class UnitigSupportClass : uint8_t {
    BackboneOnly,
    Mixed,
    ReadOnly,
};

inline bool execution_mode_runs_haplotype(ExecutionMode mode) {
    return mode == ExecutionMode::Haplotype || mode == ExecutionMode::Both;
}

inline bool execution_mode_runs_sv(ExecutionMode mode) {
    return mode == ExecutionMode::Sv || mode == ExecutionMode::Both;
}

inline UnitigSupportClass classify_unitig_support(size_t backbone_node_count,
                                                 size_t read_node_count) {
    if (backbone_node_count > 0 && read_node_count == 0) {
        return UnitigSupportClass::BackboneOnly;
    }
    if (backbone_node_count > 0) {
        return UnitigSupportClass::Mixed;
    }
    return UnitigSupportClass::ReadOnly;
}

inline const char* unitig_support_class_name(UnitigSupportClass support_class) {
    switch (support_class) {
    case UnitigSupportClass::BackboneOnly:
        return "backbone";
    case UnitigSupportClass::Mixed:
        return "mixed";
    case UnitigSupportClass::ReadOnly:
        return "read";
    }
    return "read";
}

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
    std::set<int32_t> ref_positions;

    void add_ref_pos(int32_t pos) {
        if (pos < 0) {
            return;
        }
        ref_positions.insert(pos);
        if (!is_backbone && (ref_pos < 0 || pos < ref_pos)) {
            ref_pos = pos;
        }
    }

    bool has_ref_pos_in_range(int32_t start, int32_t end) const {
        if (ref_positions.empty()) {
            return ref_pos >= start && ref_pos < end;
        }

        auto it = ref_positions.lower_bound(start);
        return it != ref_positions.end() && *it < end;
    }
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
    int32_t                  ref_id = -1;
    int32_t                  ref_start  = 0; // 0-based leftmost mapping pos
    int32_t                  ref_end    = 0; // 0-based exclusive
    int32_t                  mate_ref_id = -1;
    int32_t                  mate_ref_start = -1;
    int32_t                  template_length = 0;
    uint8_t                  mapq = 0;
    uint16_t                 flag       = 0;
    bool                     has_sa_tag = false; // supplementary alignment tag

    bool is_unmapped()      const { return flag & 0x4; }
    bool is_reverse()       const { return flag & 0x10; }
    bool mate_is_reverse()  const { return flag & 0x20; }
    bool mate_unmapped()    const { return flag & 0x8; }
    bool is_secondary()     const { return flag & 0x100; }
    bool is_supplementary() const { return flag & 0x800; }
    bool is_proper_pair()   const { return flag & 0x2; }
    bool is_read1()         const { return flag & 0x40; }
    bool mate_on_same_ref() const {
        return !is_unmapped() && !mate_unmapped() && ref_id >= 0 && ref_id == mate_ref_id;
    }
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
    int32_t               ref_pos = -1;
    std::set<int32_t>     ref_positions;
    std::vector<uint64_t> node_ids; // constituent node IDs
    size_t                backbone_node_count = 0;
    size_t                read_node_count = 0;

    UnitigSupportClass support_class() const {
        return classify_unitig_support(backbone_node_count, read_node_count);
    }
};

struct StructuralVariantCall {
    std::string chrom;
    int32_t pos = 0;
    int32_t end = 0;
    std::string id;
    std::string ref;
    std::string alt;
    std::string sv_type;
    int32_t sv_len = 0;
    std::string filter = "PASS";
    double support_score = 0.0;
    std::vector<std::string> info_fields;
};

} // namespace sharda
