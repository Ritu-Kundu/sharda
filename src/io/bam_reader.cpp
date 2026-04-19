#include "io/bam_reader.h"

#include <htslib/sam.h>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <vector>
#include <spdlog/spdlog.h>

namespace sharda {

namespace {

AlignedRead bam1_to_read(bam1_t* b, bam_hdr_t* /*hdr*/) {
    AlignedRead r;
    r.name = bam_get_qname(b);
    r.flag = b->core.flag;
    r.ref_start = b->core.pos; // 0-based

    // Sequence
    uint8_t* seq_ptr = bam_get_seq(b);
    r.seq.resize(b->core.l_qseq);
    for (int i = 0; i < b->core.l_qseq; ++i)
        r.seq[i] = seq_nt16_str[bam_seqi(seq_ptr, i)];

    // Quality
    uint8_t* qual_ptr = bam_get_qual(b);
    r.qual.resize(b->core.l_qseq);
    for (int i = 0; i < b->core.l_qseq; ++i)
        r.qual[i] = static_cast<char>(qual_ptr[i] + 33);

    // CIGAR
    uint32_t* cigar_ptr = bam_get_cigar(b);
    r.cigar.resize(b->core.n_cigar);
    int32_t ref_consumed = 0;
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        uint32_t op  = bam_cigar_op(cigar_ptr[i]);
        uint32_t len = bam_cigar_oplen(cigar_ptr[i]);
        r.cigar[i] = {static_cast<CigarOp>(op), len};
        // Track ref consumption for ref_end
        if (bam_cigar_type(op) & 2) ref_consumed += len;
    }
    r.ref_end = r.ref_start + ref_consumed;

    // SA tag
    r.has_sa_tag = (bam_aux_get(b, "SA") != nullptr);

    return r;
}

} // anonymous namespace

void iterate_read_pairs(const std::string& bam_path,
                        std::function<void(ReadPair&&)> callback) {
    samFile* fp = sam_open(bam_path.c_str(), "r");
    if (!fp) throw std::runtime_error("Cannot open BAM: " + bam_path);

    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) { sam_close(fp); throw std::runtime_error("Cannot read BAM header"); }

    bam1_t* b = bam_init1();
    AlignedRead pending;
    bool has_pending = false;

    while (sam_read1(fp, hdr, b) >= 0) {
        // Skip secondary and supplementary
        if (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        // Skip unmapped
        if (b->core.flag & BAM_FUNMAP) continue;

        AlignedRead cur = bam1_to_read(b, hdr);

        if (has_pending && pending.name == cur.name) {
            // Pair found
            ReadPair pair;
            if (pending.is_read1()) {
                pair.read1 = std::move(pending);
                pair.read2 = std::move(cur);
            } else {
                pair.read1 = std::move(cur);
                pair.read2 = std::move(pending);
            }
            has_pending = false;
            callback(std::move(pair));
        } else {
            if (has_pending) {
                spdlog::debug("Unpaired read skipped: {}", pending.name);
            }
            pending = std::move(cur);
            has_pending = true;
        }
    }
    if (has_pending) {
        spdlog::debug("Unpaired read skipped: {}", pending.name);
    }

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);
}

void create_region_bam(const std::string& src_bam,
                       const std::string& chrom,
                       int32_t start, int32_t end,
                       const std::string& dest_bam) {
    // Open source BAM
    samFile* in = sam_open(src_bam.c_str(), "r");
    if (!in)
        throw std::runtime_error("Cannot open BAM: " + src_bam);

    bam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) {
        sam_close(in);
        throw std::runtime_error("Cannot read BAM header: " + src_bam);
    }

    // Load BAM index
    hts_idx_t* idx = sam_index_load(in, src_bam.c_str());
    if (!idx) {
        bam_hdr_destroy(hdr);
        sam_close(in);
        throw std::runtime_error("Cannot load BAM index for: " + src_bam
                                 + " (is the .bai file present?)");
    }

    // Build region query string: "chrom:start+1-end" (htslib uses 1-based)
    std::string region_str = chrom + ":" + std::to_string(start + 1)
                           + "-" + std::to_string(end);

    hts_itr_t* iter = sam_itr_querys(idx, hdr, region_str.c_str());
    if (!iter) {
        hts_idx_destroy(idx);
        bam_hdr_destroy(hdr);
        sam_close(in);
        throw std::runtime_error("Cannot query region: " + region_str);
    }

    // Collect all records for the region
    std::vector<bam1_t*> records;
    bam1_t* b = bam_init1();
    while (sam_itr_next(in, iter, b) >= 0) {
        // Skip secondary, supplementary, unmapped
        if (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP))
            continue;
        bam1_t* copy = bam_dup1(b);
        records.push_back(copy);
    }
    bam_destroy1(b);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_close(in);

    spdlog::debug("Region {} extracted {} primary reads", region_str, records.size());

    // Sort by query name
    std::sort(records.begin(), records.end(), [](const bam1_t* a, const bam1_t* b) {
        return strcmp(bam_get_qname(a), bam_get_qname(b)) < 0;
    });

    // Write to destination BAM
    samFile* out = sam_open(dest_bam.c_str(), "wb");
    if (!out) {
        for (auto* r : records) bam_destroy1(r);
        bam_hdr_destroy(hdr);
        throw std::runtime_error("Cannot create BAM: " + dest_bam);
    }

    if (sam_hdr_write(out, hdr) < 0) {
        sam_close(out);
        for (auto* r : records) bam_destroy1(r);
        bam_hdr_destroy(hdr);
        throw std::runtime_error("Cannot write BAM header: " + dest_bam);
    }

    for (auto* r : records) {
        if (sam_write1(out, hdr, r) < 0) {
            sam_close(out);
            for (auto* rec : records) bam_destroy1(rec);
            bam_hdr_destroy(hdr);
            throw std::runtime_error("Cannot write BAM record: " + dest_bam);
        }
        bam_destroy1(r);
    }

    sam_close(out);
    bam_hdr_destroy(hdr);
}

} // namespace sharda
