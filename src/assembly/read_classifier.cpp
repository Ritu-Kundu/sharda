#include "assembly/read_classifier.h"

namespace sharda {

namespace {

bool cigar_has_evidence(const std::vector<CigarElement>& cigar) {
    for (const auto& c : cigar) {
        switch (c.op) {
            case CigarOp::S:  return true; // soft-clip
            case CigarOp::I:  return true; // insertion
            case CigarOp::D:  return true; // deletion
            case CigarOp::X:  return true; // mismatch
            default: break;
        }
    }
    // Also check M ops — M can hide mismatches, but we can't tell without
    // the MD tag or reference. We accept M as non-evidence here; the
    // spec says "any mismatch or indel in its cigar" which covers I/D/X.
    return false;
}

bool overlaps_any_tr(int32_t start, int32_t end, const std::vector<TandemRepeat>& trs) {
    for (const auto& tr : trs) {
        if (start < tr.end && end > tr.start) return true;
    }
    return false;
}

/// Return the TR id if the read is completely inside a single TR, else -1.
int inside_tr(int32_t start, int32_t end, const std::vector<TandemRepeat>& trs) {
    for (const auto& tr : trs) {
        if (start >= tr.start && end <= tr.end)
            return tr.id;
    }
    return -1;
}

} // anonymous namespace

ReadClassification classify_read(const AlignedRead& read,
                                 const std::vector<TandemRepeat>& trs) {
    ReadClassification cls;

    // Evidence check
    bool evidence = false;
    evidence |= cigar_has_evidence(read.cigar);
    evidence |= read.has_sa_tag;                 // split alignment
    evidence |= !read.is_proper_pair();           // improper pair
    evidence |= overlaps_any_tr(read.ref_start, read.ref_end, trs);
    cls.is_evidence = evidence;

    // ORR vs IRR: IRR if completely inside a single TR
    int tr_id = inside_tr(read.ref_start, read.ref_end, trs);
    if (tr_id >= 0) {
        cls.type  = ReadType::IRR;
        cls.tr_id = tr_id;
    } else {
        cls.type  = ReadType::ORR;
        cls.tr_id = -1;
    }

    return cls;
}

} // namespace sharda
