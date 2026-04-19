#include <gtest/gtest.h>
#include "util/kmer.h"
#include "io/fasta_reader.h"
#include "io/bed_reader.h"
#include "io/fasta_writer.h"
#include "graph/types.h"
#include "graph/dbg.h"
#include "graph/backbone.h"
#include "graph/unitig_graph.h"
#include "assembly/read_classifier.h"
#include "assembly/anchor_chain.h"
#include "assembly/graph_cleaner.h"

#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

// ── Kmer tests ──────────────────────────────────────────────────────────────

TEST(Kmer, ExtractBasic) {
    auto kmers = sharda::extract_kmers("ACGTACGT", 4);
    ASSERT_EQ(kmers.size(), 5u);
    EXPECT_EQ(kmers[0], "ACGT");
    EXPECT_EQ(kmers[4], "ACGT");
}

TEST(Kmer, TooShort) {
    auto kmers = sharda::extract_kmers("ACG", 4);
    EXPECT_TRUE(kmers.empty());
}

TEST(Kmer, ExactK) {
    auto kmers = sharda::extract_kmers("ACGT", 4);
    ASSERT_EQ(kmers.size(), 1u);
    EXPECT_EQ(kmers[0], "ACGT");
}

// ── FASTA reader/writer tests ───────────────────────────────────────────────

class TempFileTest : public ::testing::Test {
protected:
    std::string tmp_dir;
    void SetUp() override {
        tmp_dir = fs::temp_directory_path() / "sharda_test";
        fs::create_directories(tmp_dir);
    }
    void TearDown() override {
        fs::remove_all(tmp_dir);
    }
    std::string tmp_path(const std::string& name) {
        return (fs::path(tmp_dir) / name).string();
    }
};

TEST_F(TempFileTest, FastaRoundTrip) {
    std::string path = tmp_path("test.fa");
    {
        std::ofstream out(path);
        out << ">seq1\nACGTACGT\nAAAA\n";
    }
    auto [name, seq] = sharda::read_fasta(path);
    EXPECT_EQ(name, "seq1");
    EXPECT_EQ(seq, "ACGTACGTAAAA");
}

TEST_F(TempFileTest, FastaWriter) {
    std::string path = tmp_path("out.fa");
    sharda::write_fasta(path, {{"hap1", "ACGTACGT"}, {"hap2", "TTTTAAAA"}});

    std::ifstream in(path);
    std::string content((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
    EXPECT_NE(content.find(">hap1"), std::string::npos);
    EXPECT_NE(content.find("ACGTACGT"), std::string::npos);
    EXPECT_NE(content.find(">hap2"), std::string::npos);
}

// ── BED reader test ─────────────────────────────────────────────────────────

TEST_F(TempFileTest, BedReader) {
    std::string path = tmp_path("test.bed");
    {
        std::ofstream out(path);
        out << "chr1\t100\t200\n";
        out << "chr1\t500\t600\n";
    }
    auto trs = sharda::read_bed(path);
    ASSERT_EQ(trs.size(), 2u);
    EXPECT_EQ(trs[0].chrom, "chr1");
    EXPECT_EQ(trs[0].start, 100);
    EXPECT_EQ(trs[0].end, 200);
    EXPECT_EQ(trs[0].id, 0);
    EXPECT_EQ(trs[1].id, 1);
}

// ── Backbone test ───────────────────────────────────────────────────────────

TEST(Backbone, SmallReference) {
    // k=3, ref="ACGTAC" → 4 kmers: ACG, CGT, GTA, TAC
    sharda::DBG graph(3);
    std::vector<sharda::TandemRepeat> trs;
    sharda::build_backbone(graph, "ACGTAC", trs);

    EXPECT_EQ(graph.node_count(), 4u);
    EXPECT_EQ(graph.edge_count(), 3u);

    // Check backbone nodes
    EXPECT_NE(graph.backbone_node_at(0), UINT64_MAX);
    EXPECT_NE(graph.backbone_node_at(3), UINT64_MAX);
    EXPECT_EQ(graph.backbone_node_at(4), UINT64_MAX); // only 4 nodes (pos 0-3)
}

TEST(Backbone, WithTR) {
    // k=3, ref="ACGTACGT" → 6 nodes, TR at [2,5)
    sharda::DBG graph(3);
    std::vector<sharda::TandemRepeat> trs = {{"chr1", 2, 5, 0}};
    sharda::build_backbone(graph, "ACGTACGT", trs);

    // Nodes at pos 2,3,4 should have tr_id=0
    auto nid2 = graph.backbone_node_at(2);
    auto nid3 = graph.backbone_node_at(3);
    auto nid4 = graph.backbone_node_at(4);
    ASSERT_NE(nid2, UINT64_MAX);
    EXPECT_EQ(graph.node(nid2).tr_id, 0);
    EXPECT_EQ(graph.node(nid3).tr_id, 0);
    EXPECT_EQ(graph.node(nid4).tr_id, 0);

    // TR nodes list
    const auto& tr_nodes = graph.tr_nodes(0);
    EXPECT_EQ(tr_nodes.size(), 3u);
}

// ── DBG basic operations test ───────────────────────────────────────────────

TEST(DBG, ReadNodes) {
    sharda::DBG graph(3);
    auto id1 = graph.add_read_node("ACG");
    auto id2 = graph.add_read_node("CGT");
    auto id3 = graph.add_read_node("ACG"); // duplicate, should return id1

    EXPECT_EQ(id3, id1);
    EXPECT_NE(id1, id2);
    EXPECT_EQ(graph.node(id1).kmer, "ACG");
    EXPECT_FALSE(graph.node(id1).is_backbone);
}

TEST(DBG, EdgeIncrement) {
    sharda::DBG graph(3);
    auto a = graph.add_read_node("ACG");
    auto b = graph.add_read_node("CGT");
    graph.add_edge(a, b);
    graph.add_edge(a, b);

    ASSERT_EQ(graph.edge_count(), 1u);
    EXPECT_EQ(graph.edges()[0].weight, 2u);
}

TEST(DBG, HaplotypeEdges) {
    sharda::DBG graph(3);
    auto a = graph.add_read_node("ACG");
    auto b = graph.add_read_node("CGT");
    graph.add_haplotype_edge(a, b);

    ASSERT_EQ(graph.haplotype_edges().size(), 1u);
    EXPECT_EQ(graph.haplotype_edges()[0].from_node, a);
    EXPECT_EQ(graph.haplotype_edges()[0].to_node, b);
}

// ── Read classifier test ────────────────────────────────────────────────────

TEST(ReadClassifier, EvidenceSoftClip) {
    sharda::AlignedRead read;
    read.seq = "ACGTACGT";
    read.cigar = {{sharda::CigarOp::S, 2}, {sharda::CigarOp::M, 6}};
    read.ref_start = 10;
    read.ref_end   = 16;
    read.flag       = 0x3; // proper pair

    std::vector<sharda::TandemRepeat> trs;
    auto cls = sharda::classify_read(read, trs);
    EXPECT_TRUE(cls.is_evidence); // soft-clip
    EXPECT_EQ(cls.type, sharda::ReadType::ORR);
}

TEST(ReadClassifier, IRRInsideTR) {
    sharda::AlignedRead read;
    read.seq = "ACGTACGT";
    read.cigar = {{sharda::CigarOp::M, 8}};
    read.ref_start = 100;
    read.ref_end   = 108;
    read.flag       = 0x3; // proper pair

    std::vector<sharda::TandemRepeat> trs = {{"chr1", 50, 200, 0}};
    auto cls = sharda::classify_read(read, trs);
    EXPECT_TRUE(cls.is_evidence); // overlaps TR
    EXPECT_EQ(cls.type, sharda::ReadType::IRR);
    EXPECT_EQ(cls.tr_id, 0);
}

TEST(ReadClassifier, ORRNoEvidence) {
    sharda::AlignedRead read;
    read.seq = "ACGTACGT";
    read.cigar = {{sharda::CigarOp::M, 8}};
    read.ref_start = 10;
    read.ref_end   = 18;
    read.flag       = 0x3; // proper pair

    std::vector<sharda::TandemRepeat> trs = {{"chr1", 500, 600, 0}};
    auto cls = sharda::classify_read(read, trs);
    EXPECT_FALSE(cls.is_evidence);
    EXPECT_EQ(cls.type, sharda::ReadType::ORR);
}

// ── Unitig graph test ───────────────────────────────────────────────────────

TEST(UnitigGraph, LinearCollapse) {
    // Build a linear graph: A -> B -> C -> D
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTAC", {}); // 4 nodes, 3 edges

    sharda::UnitigGraph ug;
    bool ok = ug.build(graph);
    EXPECT_TRUE(ok);
    // Should collapse into 1 unitig
    EXPECT_EQ(ug.unitig_count(), 1u);
    EXPECT_EQ(ug.edges().size(), 0u);
}

TEST(UnitigGraph, BranchPreserved) {
    // Build a graph with a branch
    sharda::DBG graph(3);
    auto a = graph.add_backbone_node("ACG", 0, -1);
    auto b = graph.add_backbone_node("CGT", 1, -1);
    auto c = graph.add_backbone_node("GTA", 2, -1);
    auto d = graph.add_read_node("GTT"); // branch at position 2

    graph.add_edge(a, b);
    graph.add_edge(b, c);
    graph.add_edge(b, d);

    sharda::UnitigGraph ug;
    bool ok = ug.build(graph);
    EXPECT_TRUE(ok);
    EXPECT_GT(ug.unitig_count(), 1u); // should not collapse everything
}

// ── Target regions BED reader test ──────────────────────────────────────────

TEST_F(TempFileTest, TargetRegionsReader) {
    std::string path = tmp_path("targets.bed");
    {
        std::ofstream out(path);
        out << "chr1\t10000\t20000\n";
        out << "chr2\t50000\t60000\n";
    }
    auto regions = sharda::read_target_regions(path);
    ASSERT_EQ(regions.size(), 2u);
    EXPECT_EQ(regions[0].chrom, "chr1");
    EXPECT_EQ(regions[0].start, 10000);
    EXPECT_EQ(regions[0].end, 20000);
    EXPECT_EQ(regions[1].chrom, "chr2");
}

// ── TR filtering for region test ────────────────────────────────────────────

TEST(TRFilter, FiltersAndAdjusts) {
    std::vector<sharda::TandemRepeat> all_trs = {
        {"chr1", 12000, 13000, 0},
        {"chr1", 50000, 51000, 1},  // outside target region
        {"chr2", 12000, 13000, 2},  // wrong chrom
    };
    sharda::TargetRegion region{"chr1", 10000, 20000};
    int32_t padding = 1000;

    auto local = sharda::filter_trs_for_region(all_trs, region, padding);
    ASSERT_EQ(local.size(), 1u);
    // ext_start = 10000 - 1000 = 9000, so local start = 12000 - 9000 = 3000
    EXPECT_EQ(local[0].start, 3000);
    EXPECT_EQ(local[0].end, 4000);
    EXPECT_EQ(local[0].id, 0);
}

TEST(TRFilter, ClampsToBoundary) {
    std::vector<sharda::TandemRepeat> all_trs = {
        {"chr1", 500, 2000, 0},  // starts before ext_start=0
    };
    sharda::TargetRegion region{"chr1", 1000, 5000};
    int32_t padding = 1000;
    // ext_start = max(0, 1000-1000) = 0

    auto local = sharda::filter_trs_for_region(all_trs, region, padding);
    ASSERT_EQ(local.size(), 1u);
    EXPECT_EQ(local[0].start, 500);  // max(500, 0) - 0 = 500
    EXPECT_EQ(local[0].end, 2000);   // min(2000, 6000) - 0 = 2000
}
