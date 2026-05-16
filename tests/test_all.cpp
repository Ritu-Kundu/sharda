#include "io/vcf_writer.h"
#include <gtest/gtest.h>
#include "util/kmer.h"
#include "io/fasta_reader.h"
#include "io/bed_reader.h"
#include "io/debug_artifacts.h"
#include "io/fasta_writer.h"
#include "io/gfa_writer.h"
#include "util/debug_query.h"
#include "graph/types.h"
#include "graph/dbg.h"
#include "graph/backbone.h"
#include "graph/unitig_graph.h"
#include "assembly/read_classifier.h"
#include "assembly/anchor_chain.h"
#include "assembly/read_adder.h"
#include "assembly/graph_cleaner.h"
#include "assembly/flow_decomp.h"
#include "assembly/sv_caller.h"

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

namespace {

void add_edge_copies(sharda::DBG& graph, uint64_t from, uint64_t to, uint32_t copies) {
    for (uint32_t i = 0; i < copies; ++i) {
        graph.add_edge(from, to);
    }
}

bool has_edge(const sharda::DBG& graph, uint64_t from, uint64_t to) {
    for (uint64_t edge_index : graph.out_edges(from)) {
        if (graph.edges()[edge_index].to == to) {
            return true;
        }
    }
    return false;
}

std::string read_text_file(const std::string& path) {
    std::ifstream in(path);
    return std::string((std::istreambuf_iterator<char>(in)),
                       std::istreambuf_iterator<char>());
}

} // namespace

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

TEST(GraphCleaner, RemovesUnderSupportedTipEvenWithHaplotypeEdges) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTACGT", {});

    auto b2 = graph.backbone_node_at(2);
    auto b3 = graph.backbone_node_at(3);
    auto b4 = graph.backbone_node_at(4);
    auto b5 = graph.backbone_node_at(5);
    ASSERT_NE(b2, UINT64_MAX);
    ASSERT_NE(b3, UINT64_MAX);
    ASSERT_NE(b4, UINT64_MAX);
    ASSERT_NE(b5, UINT64_MAX);

    add_edge_copies(graph, b2, b3, 99);
    add_edge_copies(graph, b3, b4, 99);
    add_edge_copies(graph, b4, b5, 99);

    auto t1 = graph.add_read_node("TTT");
    auto t2 = graph.add_read_node("TTA");
    graph.add_node_ref_pos(t1, 2);
    graph.add_node_ref_pos(t2, 3);
    graph.add_edge(t1, t2);
    graph.add_edge(t2, b4);
    graph.add_haplotype_edge(t1, t2);

    sharda::clean_graph(graph, 150);

    EXPECT_TRUE(graph.is_node_removed(t1));
    EXPECT_TRUE(graph.is_node_removed(t2));
    EXPECT_TRUE(graph.haplotype_edges().empty());
}

TEST(GraphCleaner, KeepsSupportedShortTip) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTACGT", {});

    auto b2 = graph.backbone_node_at(2);
    auto b3 = graph.backbone_node_at(3);
    auto b4 = graph.backbone_node_at(4);
    auto b5 = graph.backbone_node_at(5);
    ASSERT_NE(b2, UINT64_MAX);
    ASSERT_NE(b3, UINT64_MAX);
    ASSERT_NE(b4, UINT64_MAX);
    ASSERT_NE(b5, UINT64_MAX);

    add_edge_copies(graph, b2, b3, 99);
    add_edge_copies(graph, b3, b4, 99);
    add_edge_copies(graph, b4, b5, 99);

    auto t1 = graph.add_read_node("GGG");
    auto t2 = graph.add_read_node("GGA");
    graph.add_node_ref_pos(t1, 2);
    graph.add_node_ref_pos(t2, 3);
    add_edge_copies(graph, t1, t2, 10);
    add_edge_copies(graph, t2, b4, 10);
    graph.add_haplotype_edge(t1, t2);

    sharda::clean_graph(graph, 150);

    EXPECT_FALSE(graph.is_node_removed(t1));
    EXPECT_FALSE(graph.is_node_removed(t2));
    EXPECT_TRUE(has_edge(graph, t1, t2));
    EXPECT_TRUE(has_edge(graph, t2, b4));
}

TEST(GraphCleaner, PrunesLowWeightNonBackboneEdgesUsingImpliedCoordinates) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTACGT", {});

    auto b0 = graph.backbone_node_at(0);
    auto b1 = graph.backbone_node_at(1);
    auto b2 = graph.backbone_node_at(2);
    auto b3 = graph.backbone_node_at(3);
    ASSERT_NE(b0, UINT64_MAX);
    ASSERT_NE(b1, UINT64_MAX);
    ASSERT_NE(b2, UINT64_MAX);
    ASSERT_NE(b3, UINT64_MAX);

    add_edge_copies(graph, b0, b1, 99);
    add_edge_copies(graph, b1, b2, 99);
    add_edge_copies(graph, b2, b3, 99);

    auto r1 = graph.add_read_node("TTT");
    auto r2 = graph.add_read_node("TTC");
    graph.add_node_ref_pos(r1, 1);
    graph.add_node_ref_pos(r2, 2);

    add_edge_copies(graph, b0, r1, 10);
    graph.add_edge(r1, r2);
    add_edge_copies(graph, r2, b3, 10);

    sharda::clean_graph(graph, 150);

    EXPECT_FALSE(has_edge(graph, r1, r2));
    EXPECT_TRUE(has_edge(graph, b0, r1));
    EXPECT_TRUE(has_edge(graph, r2, b3));
}

TEST(FlowDecomposition, UsesAnchoredBoundaryUnitigsWhenTopologyIsAmbiguous) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTACGT", {});

    auto extra1 = graph.add_read_node("TTT");
    auto extra2 = graph.add_read_node("TTA");
    graph.add_edge(extra1, extra2);

    sharda::UnitigGraph ug;
    ASSERT_TRUE(ug.build(graph));
    ASSERT_EQ(ug.unitig_count(), 2u);

    auto unanchored_paths = sharda::flow_decomposition(ug, 2);
    EXPECT_TRUE(unanchored_paths.empty());

    sharda::FlowBoundaryAnchors anchors;
    anchors.start_node_id = graph.backbone_node_at(0);
    anchors.end_node_id = graph.backbone_node_at(5);

    auto anchored_paths = sharda::flow_decomposition(ug, 2, anchors);
    ASSERT_EQ(anchored_paths.size(), 1u);
    ASSERT_EQ(anchored_paths[0].unitig_ids.size(), 1u);

    const uint64_t anchored_unitig = ug.node_to_unitig(anchors.start_node_id);
    EXPECT_EQ(anchored_unitig, ug.node_to_unitig(anchors.end_node_id));
    EXPECT_EQ(anchored_paths[0].unitig_ids.front(), anchored_unitig);
    EXPECT_EQ(anchored_paths[0].sequence, ug.unitig(anchored_unitig).sequence);
}

TEST(GraphCleaner, DoesNotPopBubblesDuringCleaning) {
    sharda::DBG graph(3);
    auto source = graph.add_backbone_node("AAA", 0, -1);
    auto sink = graph.add_backbone_node("CCC", 3, -1);
    auto strong = graph.add_read_node("AAT");
    auto weak = graph.add_read_node("AAC");

    graph.add_node_ref_pos(strong, 1);
    graph.add_node_ref_pos(weak, 1);

    add_edge_copies(graph, source, strong, 19);
    add_edge_copies(graph, strong, sink, 19);
    add_edge_copies(graph, source, weak, 3);
    add_edge_copies(graph, weak, sink, 3);

    sharda::clean_graph(graph, 150);

    EXPECT_FALSE(graph.is_node_removed(strong));
    EXPECT_FALSE(graph.is_node_removed(weak));
    EXPECT_TRUE(has_edge(graph, source, weak));
    EXPECT_TRUE(has_edge(graph, weak, sink));
}

TEST(GraphCleaner, RemovesTipUsingRegionalMeanDepthFloor) {
    sharda::DBG graph(3);
    auto b0 = graph.add_backbone_node("AAA", 0, -1);
    auto b1 = graph.add_backbone_node("AAT", 1, -1);
    auto b2 = graph.add_backbone_node("ATC", 2, -1);
    graph.node_mut(b0).depth = 8;
    graph.node_mut(b1).depth = 8;
    graph.node_mut(b2).depth = 8;

    add_edge_copies(graph, b0, b1, 4);
    add_edge_copies(graph, b1, b2, 4);

    auto t0 = graph.add_read_node("CCA");
    auto t1 = graph.add_read_node("CCC");
    graph.add_node_ref_pos(t0, 0);
    graph.add_node_ref_pos(t1, 1);
    graph.add_edge(b0, t0);
    graph.add_edge(t0, t1);

    sharda::clean_graph(graph, 150);

    EXPECT_TRUE(graph.is_node_removed(t0));
    EXPECT_TRUE(graph.is_node_removed(t1));
}

TEST(GraphCleaner, PrunesEdgesUsingRegionalMeanDepthFloor) {
    sharda::DBG graph(3);
    auto b0 = graph.add_backbone_node("AAA", 0, -1);
    auto b1 = graph.add_backbone_node("AAT", 1, -1);
    auto b2 = graph.add_backbone_node("ATC", 2, -1);
    graph.node_mut(b0).depth = 8;
    graph.node_mut(b1).depth = 8;
    graph.node_mut(b2).depth = 8;

    add_edge_copies(graph, b0, b1, 4);
    add_edge_copies(graph, b1, b2, 4);

    auto r0 = graph.add_read_node("CCA");
    auto r1 = graph.add_read_node("CCC");
    graph.add_node_ref_pos(r0, 0);
    graph.add_node_ref_pos(r1, 1);
    add_edge_copies(graph, b0, r0, 2);
    graph.add_edge(r0, r1);
    add_edge_copies(graph, r1, b2, 2);

    sharda::clean_graph(graph, 150);

    EXPECT_FALSE(has_edge(graph, r0, r1));
    EXPECT_TRUE(has_edge(graph, b0, r0));
    EXPECT_TRUE(has_edge(graph, r1, b2));
}

TEST(GraphCleaner, PrunesWeakInternalAlternateBranch) {
    sharda::DBG graph(3);
    auto b0 = graph.add_backbone_node("AAA", 0, -1);
    auto b1 = graph.add_backbone_node("AAT", 1, -1);
    auto b2 = graph.add_backbone_node("ATC", 2, -1);
    auto b3 = graph.add_backbone_node("TCG", 3, -1);
    auto b4 = graph.add_backbone_node("CGT", 4, -1);

    graph.node_mut(b0).depth = 4;
    graph.node_mut(b1).depth = 4;
    graph.node_mut(b2).depth = 4;
    graph.node_mut(b3).depth = 4;
    graph.node_mut(b4).depth = 4;

    add_edge_copies(graph, b0, b1, 4);
    add_edge_copies(graph, b1, b2, 4);
    add_edge_copies(graph, b2, b3, 4);
    add_edge_copies(graph, b3, b4, 4);

    auto r0 = graph.add_read_node("CCA");
    auto r1 = graph.add_read_node("CCC");
    graph.add_node_ref_pos(r0, 1);
    graph.add_node_ref_pos(r1, 2);

    graph.add_edge(b1, r0);
    graph.add_edge(r0, r1);
    graph.add_edge(r1, b3);

    sharda::clean_graph(graph, 150);

    EXPECT_TRUE(graph.is_node_removed(r0));
    EXPECT_TRUE(graph.is_node_removed(r1));
    EXPECT_FALSE(has_edge(graph, b1, r0));
    EXPECT_FALSE(has_edge(graph, r0, r1));
    EXPECT_FALSE(has_edge(graph, r1, b3));
    EXPECT_TRUE(has_edge(graph, b1, b2));
    EXPECT_TRUE(has_edge(graph, b2, b3));
}

TEST(GraphCleaner, PrunesWeakInternalAlternateComponentWithSplit) {
    sharda::DBG graph(3);
    auto b0 = graph.add_backbone_node("AAA", 0, -1);
    auto b1 = graph.add_backbone_node("AAT", 1, -1);
    auto b2 = graph.add_backbone_node("ATC", 2, -1);
    auto b3 = graph.add_backbone_node("TCG", 3, -1);
    auto b4 = graph.add_backbone_node("CGT", 4, -1);

    graph.node_mut(b0).depth = 4;
    graph.node_mut(b1).depth = 4;
    graph.node_mut(b2).depth = 4;
    graph.node_mut(b3).depth = 4;
    graph.node_mut(b4).depth = 4;

    add_edge_copies(graph, b0, b1, 4);
    add_edge_copies(graph, b1, b2, 4);
    add_edge_copies(graph, b2, b3, 4);
    add_edge_copies(graph, b3, b4, 4);

    auto r0 = graph.add_read_node("CCA");
    auto r1 = graph.add_read_node("CCC");
    auto r2 = graph.add_read_node("CCG");
    auto r3 = graph.add_read_node("CCT");
    graph.add_node_ref_pos(r0, 1);
    graph.add_node_ref_pos(r1, 2);
    graph.add_node_ref_pos(r2, 3);
    graph.add_node_ref_pos(r3, 3);

    graph.add_edge(b1, r0);
    graph.add_edge(r0, r1);
    graph.add_edge(r1, r2);
    graph.add_edge(r1, r3);
    graph.add_edge(r2, b4);
    graph.add_edge(r3, b4);

    sharda::clean_graph(graph, 150);

    EXPECT_TRUE(graph.is_node_removed(r0));
    EXPECT_TRUE(graph.is_node_removed(r1));
    EXPECT_TRUE(graph.is_node_removed(r2));
    EXPECT_TRUE(graph.is_node_removed(r3));
    EXPECT_FALSE(has_edge(graph, b1, r0));
    EXPECT_FALSE(has_edge(graph, r0, r1));
    EXPECT_FALSE(has_edge(graph, r1, r2));
    EXPECT_FALSE(has_edge(graph, r1, r3));
    EXPECT_FALSE(has_edge(graph, r2, b4));
    EXPECT_FALSE(has_edge(graph, r3, b4));
    EXPECT_TRUE(has_edge(graph, b1, b2));
    EXPECT_TRUE(has_edge(graph, b2, b3));
    EXPECT_TRUE(has_edge(graph, b3, b4));
}

TEST(GraphCleaner, PreservesBackboneEdgesInSvMode) {
    sharda::DBG default_graph(3);
    auto default_b0 = default_graph.add_backbone_node("AAA", 0, -1);
    auto default_b1 = default_graph.add_backbone_node("AAT", 1, -1);
    auto default_b2 = default_graph.add_backbone_node("ATC", 2, -1);
    auto default_b3 = default_graph.add_backbone_node("TCG", 3, -1);

    add_edge_copies(default_graph, default_b0, default_b1, 100);
    default_graph.add_edge(default_b1, default_b2);
    add_edge_copies(default_graph, default_b2, default_b3, 100);

    sharda::clean_graph(default_graph, 150);
    EXPECT_FALSE(has_edge(default_graph, default_b1, default_b2));

    sharda::DBG sv_graph(3);
    auto sv_b0 = sv_graph.add_backbone_node("AAA", 0, -1);
    auto sv_b1 = sv_graph.add_backbone_node("AAT", 1, -1);
    auto sv_b2 = sv_graph.add_backbone_node("ATC", 2, -1);
    auto sv_b3 = sv_graph.add_backbone_node("TCG", 3, -1);

    add_edge_copies(sv_graph, sv_b0, sv_b1, 100);
    sv_graph.add_edge(sv_b1, sv_b2);
    add_edge_copies(sv_graph, sv_b2, sv_b3, 100);

    sharda::GraphCleaningOptions options;
    options.preserve_backbone_edges = true;
    sharda::clean_graph(sv_graph, 150, options);
    EXPECT_TRUE(has_edge(sv_graph, sv_b1, sv_b2));
}

TEST_F(TempFileTest, DebugArtifactsUseStableNodeIds) {
    sharda::DBG graph_a(3);
    auto a0 = graph_a.add_backbone_node("AAA", 0, -1);
    auto a1 = graph_a.add_backbone_node("AAT", 1, -1);
    auto ax = graph_a.add_read_node("CCA");
    auto ay = graph_a.add_read_node("CCG");
    graph_a.add_node_ref_pos(ax, 1);
    graph_a.add_node_ref_pos(ay, 2);
    graph_a.add_edge(a0, ax);
    graph_a.add_edge(ax, ay);
    graph_a.add_edge(ay, a1);

    sharda::DBG graph_b(3);
    auto b0 = graph_b.add_backbone_node("AAA", 0, -1);
    auto b1 = graph_b.add_backbone_node("AAT", 1, -1);
    auto by = graph_b.add_read_node("CCG");
    auto bx = graph_b.add_read_node("CCA");
    graph_b.add_node_ref_pos(bx, 1);
    graph_b.add_node_ref_pos(by, 2);
    graph_b.add_edge(b0, bx);
    graph_b.add_edge(bx, by);
    graph_b.add_edge(by, b1);

    const std::string gfa_a = tmp_path("graph_a.gfa");
    const std::string gfa_b = tmp_path("graph_b.gfa");
    const std::string json_a = tmp_path("graph_a.json");
    const std::string json_b = tmp_path("graph_b.json");
    sharda::write_gfa(gfa_a, graph_a);
    sharda::write_gfa(gfa_b, graph_b);
    sharda::write_dbg_json(json_a, graph_a);
    sharda::write_dbg_json(json_b, graph_b);

    EXPECT_EQ(read_text_file(gfa_a), read_text_file(gfa_b));
    EXPECT_EQ(read_text_file(json_a), read_text_file(json_b));
}

TEST_F(TempFileTest, SvUnitigArtifactsAnnotateSupportClass) {
    sharda::DBG graph(3);

    auto mixed_backbone = graph.add_backbone_node("AAA", 0, -1);
    auto mixed_read = graph.add_read_node("AAT");
    graph.add_node_ref_pos(mixed_read, 1);
    graph.add_edge(mixed_backbone, mixed_read);

    auto pure_backbone_a = graph.add_backbone_node("ACC", 3, -1);
    auto pure_backbone_b = graph.add_backbone_node("CCG", 4, -1);
    graph.add_edge(pure_backbone_a, pure_backbone_b);

    auto read_only_a = graph.add_read_node("TTT");
    auto read_only_b = graph.add_read_node("TTC");
    graph.add_node_ref_pos(read_only_a, 6);
    graph.add_node_ref_pos(read_only_b, 7);
    graph.add_edge(read_only_a, read_only_b);

    sharda::UnitigGraph ug;
    ASSERT_TRUE(ug.build(graph));

    std::string gfa_path = tmp_path("unitig.sv.gfa");
    std::string json_path = tmp_path("unitig.sv.json");
    sharda::write_sv_unitig_gfa(gfa_path, ug);
    sharda::write_sv_unitig_json(json_path, ug);

    const std::string gfa = read_text_file(gfa_path);
    const std::string json = read_text_file(json_path);

    EXPECT_NE(gfa.find("SC:Z:backbone"), std::string::npos);
    EXPECT_NE(gfa.find("SC:Z:mixed"), std::string::npos);
    EXPECT_NE(gfa.find("SC:Z:read"), std::string::npos);
    EXPECT_NE(gfa.find("CL:z:#3B7A57"), std::string::npos);
    EXPECT_NE(json.find("\"graph_kind\": \"unitig_sv\""), std::string::npos);
    EXPECT_NE(json.find("\"support_class\": \"backbone\""), std::string::npos);
    EXPECT_NE(json.find("\"support_class\": \"mixed\""), std::string::npos);
    EXPECT_NE(json.find("\"support_class\": \"read\""), std::string::npos);
}

TEST(SvCaller, CallsMultipleParallelInsertionPathsAgainstBackbonePath) {
    sharda::DBG graph(3);

    auto source = graph.add_backbone_node("AAA", 0, -1);
    auto ref1 = graph.add_backbone_node("AAC", 1, -1);
    auto ref2 = graph.add_backbone_node("ACC", 2, -1);
    auto sink = graph.add_backbone_node("CCC", 3, -1);

    graph.add_edge(source, ref1);
    graph.add_edge(ref1, ref2);
    graph.add_edge(ref2, sink);

    auto alt1_a = graph.add_read_node("AAG");
    auto alt1_b = graph.add_read_node("AGC");
    auto alt1_c = graph.add_read_node("GCC");
    graph.add_node_ref_pos(alt1_a, 1);
    graph.add_node_ref_pos(alt1_b, 2);
    graph.add_node_ref_pos(alt1_c, 3);
    graph.add_edge(source, alt1_a);
    graph.add_edge(alt1_a, alt1_b);
    graph.add_edge(alt1_b, alt1_c);
    graph.add_edge(alt1_c, sink);

    auto alt2_a = graph.add_read_node("AAT");
    auto alt2_b = graph.add_read_node("ATC");
    auto alt2_c = graph.add_read_node("TCC");
    graph.add_node_ref_pos(alt2_a, 1);
    graph.add_node_ref_pos(alt2_b, 2);
    graph.add_node_ref_pos(alt2_c, 3);
    graph.add_edge(source, alt2_a);
    graph.add_edge(alt2_a, alt2_b);
    graph.add_edge(alt2_b, alt2_c);
    graph.add_edge(alt2_c, sink);

    sharda::UnitigGraph ug;
    ASSERT_TRUE(ug.build(graph));

    auto calls = sharda::call_structural_variants(ug, "chrTest", 100);
    ASSERT_EQ(calls.size(), 2u);

    std::sort(calls.begin(), calls.end(), [](const sharda::StructuralVariantCall& lhs,
                                             const sharda::StructuralVariantCall& rhs) {
        return lhs.alt < rhs.alt;
    });

    EXPECT_EQ(calls[0].chrom, "chrTest");
    EXPECT_EQ(calls[0].sv_type, "INS");
    EXPECT_EQ(calls[0].ref, "A");
    EXPECT_EQ(calls[0].alt, "AG");
    EXPECT_EQ(calls[0].pos, 102);
    EXPECT_EQ(calls[0].end, 103);
    EXPECT_EQ(calls[0].sv_len, 1);
    EXPECT_NE(std::find(calls[0].info_fields.begin(), calls[0].info_fields.end(), "SRC_REF_POS=101"),
              calls[0].info_fields.end());
    EXPECT_NE(std::find(calls[0].info_fields.begin(), calls[0].info_fields.end(), "SNK_REF_POS=104"),
              calls[0].info_fields.end());

    EXPECT_EQ(calls[1].sv_type, "INS");
    EXPECT_EQ(calls[1].ref, "A");
    EXPECT_EQ(calls[1].alt, "AT");
    EXPECT_EQ(calls[1].pos, 102);
    EXPECT_EQ(calls[1].end, 103);
    EXPECT_EQ(calls[1].sv_len, 1);
}

TEST(SvCaller, CallsDeletionPathAgainstBackbonePath) {
    sharda::DBG graph(3);

    auto source = graph.add_backbone_node("AAA", 0, -1);
    auto ref1 = graph.add_backbone_node("AAC", 1, -1);
    auto ref2 = graph.add_backbone_node("ACC", 2, -1);
    auto sink = graph.add_backbone_node("CCC", 3, -1);

    graph.add_edge(source, ref1);
    graph.add_edge(ref1, ref2);
    graph.add_edge(ref2, sink);

    auto alt = graph.add_read_node("ACC");
    graph.add_node_ref_pos(alt, 2);
    graph.add_edge(source, alt);
    graph.add_edge(alt, sink);

    sharda::UnitigGraph ug;
    ASSERT_TRUE(ug.build(graph));

    auto calls = sharda::call_structural_variants(ug, "chrDel", 200);
    ASSERT_EQ(calls.size(), 1u);

    EXPECT_EQ(calls[0].chrom, "chrDel");
    EXPECT_EQ(calls[0].sv_type, "DEL");
    EXPECT_EQ(calls[0].ref, "CC");
    EXPECT_EQ(calls[0].alt, "C");
    EXPECT_EQ(calls[0].pos, 204);
    EXPECT_EQ(calls[0].end, 206);
    EXPECT_EQ(calls[0].sv_len, -1);
    EXPECT_NE(std::find(calls[0].info_fields.begin(), calls[0].info_fields.end(), "SRC_REF_POS=201"),
              calls[0].info_fields.end());
    EXPECT_NE(std::find(calls[0].info_fields.begin(), calls[0].info_fields.end(), "SNK_REF_POS=204"),
              calls[0].info_fields.end());
}

TEST(SvCaller, UsesNearestBackboneRejoinForCanonicalInterval) {
    sharda::DBG graph(3);

    auto source = graph.add_backbone_node("AAA", 0, -1);
    auto ref1 = graph.add_backbone_node("AAC", 1, -1);
    auto ref2 = graph.add_backbone_node("ACC", 2, -1);
    auto sink = graph.add_backbone_node("CCC", 3, -1);
    auto tail = graph.add_backbone_node("CCG", 4, -1);

    graph.add_edge(source, ref1);
    graph.add_edge(ref1, ref2);
    graph.add_edge(ref2, sink);
    graph.add_edge(sink, tail);

    auto alt = graph.add_read_node("ACC");
    graph.add_node_ref_pos(alt, 2);
    graph.add_edge(source, alt);
    graph.add_edge(alt, sink);

    sharda::UnitigGraph ug;
    ASSERT_TRUE(ug.build(graph));

    auto calls = sharda::call_structural_variants(ug, "chrNearest", 200);
    ASSERT_EQ(calls.size(), 1u);

    EXPECT_EQ(calls[0].sv_type, "DEL");
    EXPECT_EQ(calls[0].pos, 204);
    EXPECT_EQ(calls[0].end, 206);
    EXPECT_NE(std::find(calls[0].info_fields.begin(), calls[0].info_fields.end(), "SRC_REF_POS=201"),
              calls[0].info_fields.end());
    EXPECT_NE(std::find(calls[0].info_fields.begin(), calls[0].info_fields.end(), "SNK_REF_POS=204"),
              calls[0].info_fields.end());
    EXPECT_EQ(std::find(calls[0].info_fields.begin(), calls[0].info_fields.end(), "SNK_REF_POS=205"),
              calls[0].info_fields.end());
}

TEST(SvCaller, CollapsesEquivalentCanonicalCallsBySupport) {
    sharda::StructuralVariantCall weaker;
    weaker.chrom = "chrCollapse";
    weaker.pos = 451;
    weaker.end = 453;
    weaker.ref = "TT";
    weaker.alt = "T";
    weaker.sv_type = "DEL";
    weaker.sv_len = -1;
    weaker.support_score = 3.0;
    weaker.info_fields = {"SRC_UID=0", "SNK_UID=4", "SRC_REF_POS=288", "SNK_REF_POS=453"};

    sharda::StructuralVariantCall stronger = weaker;
    stronger.support_score = 9.0;
    stronger.info_fields = {"SRC_UID=2", "SNK_UID=4", "SRC_REF_POS=288", "SNK_REF_POS=453"};

    auto collapsed = sharda::collapse_structural_variant_calls({weaker, stronger});
    ASSERT_EQ(collapsed.size(), 1u);
    EXPECT_EQ(collapsed[0].id, "sv1");
    EXPECT_EQ(collapsed[0].support_score, 9.0);
    EXPECT_NE(std::find(collapsed[0].info_fields.begin(), collapsed[0].info_fields.end(), "SRC_UID=2"),
              collapsed[0].info_fields.end());
    EXPECT_EQ(std::find(collapsed[0].info_fields.begin(), collapsed[0].info_fields.end(), "SRC_UID=0"),
              collapsed[0].info_fields.end());
}

TEST(SvCaller, CollapsesOverlappingCanonicalCallsWithSharedBoundary) {
    sharda::StructuralVariantCall earlier;
    earlier.chrom = "chrOverlap";
    earlier.pos = 451;
    earlier.end = 453;
    earlier.ref = "TT";
    earlier.alt = "T";
    earlier.sv_type = "DEL";
    earlier.sv_len = -1;
    earlier.support_score = 4.0;
    earlier.info_fields = {"SRC_UID=2", "SNK_UID=4", "SRC_REF_POS=288", "SNK_REF_POS=453"};

    sharda::StructuralVariantCall later;
    later.chrom = "chrOverlap";
    later.pos = 452;
    later.end = 454;
    later.ref = "TG";
    later.alt = "T";
    later.sv_type = "DEL";
    later.sv_len = -1;
    later.support_score = 6.5;
    later.info_fields = {"SRC_UID=3", "SNK_UID=4", "SRC_REF_POS=409", "SNK_REF_POS=453"};

    auto collapsed = sharda::collapse_structural_variant_calls({earlier, later});
    ASSERT_EQ(collapsed.size(), 1u);
    EXPECT_EQ(collapsed[0].id, "sv1");
    EXPECT_EQ(collapsed[0].pos, later.pos);
    EXPECT_EQ(collapsed[0].end, later.end);
    EXPECT_EQ(collapsed[0].support_score, later.support_score);
    EXPECT_NE(std::find(collapsed[0].info_fields.begin(), collapsed[0].info_fields.end(), "SNK_REF_POS=453"),
              collapsed[0].info_fields.end());
}

TEST_F(TempFileTest, VcfWriterIncludesSupportInfo) {
    sharda::StructuralVariantCall call;
    call.chrom = "chrSupport";
    call.pos = 10;
    call.end = 12;
    call.id = "sv1";
    call.ref = "TT";
    call.alt = "T";
    call.sv_type = "DEL";
    call.sv_len = -1;
    call.support_score = 4.25;
    call.info_fields = {"SRC_UID=1", "SNK_UID=2", "SRC_REF_POS=11", "SNK_REF_POS=13"};

    const std::string path = tmp_path("calls.vcf");
    sharda::write_vcf(path, {call}, "sharda-test");

    const std::string vcf = read_text_file(path);
    EXPECT_NE(vcf.find("##INFO=<ID=SUPPORT,Number=1,Type=Float"), std::string::npos);
    EXPECT_NE(vcf.find("SUPPORT=4.25"), std::string::npos);
}

TEST(ReadAdder, OrrChoosesClosestMatchingBackboneNodes) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTTACGTA", {});

    sharda::ReadPair pair;
    pair.read1.name = "orr_closest";
    pair.read1.seq = "ACGTA";
    pair.read1.cigar = {{sharda::CigarOp::M, 5}};
    pair.read1.ref_start = 4;
    pair.read1.ref_end = 9;
    pair.read1.flag = 0x43;

    pair.read2.name = "mate_unused";
    pair.read2.seq = "AC";
    pair.read2.cigar = {{sharda::CigarOp::M, 2}};
    pair.read2.ref_start = 0;
    pair.read2.ref_end = 2;
    pair.read2.flag = 0x83;

    std::vector<std::string> traced_reads = {"orr_closest"};
    std::vector<sharda::ReadTraceRecord> traces;
    sharda::ReadTraceSink trace_sink{&traced_reads, &traces};

    sharda::add_read_pair(pair, graph, {}, trace_sink);

    ASSERT_EQ(traces.size(), 1u);
    ASSERT_EQ(traces[0].raw_nodes.size(), 3u);
    EXPECT_EQ(traces[0].raw_nodes[0].node_id, graph.backbone_node_at(5));
    EXPECT_EQ(traces[0].raw_nodes[1].node_id, graph.backbone_node_at(6));
    EXPECT_EQ(traces[0].raw_nodes[2].node_id, graph.backbone_node_at(7));
}

TEST(ReadAdder, OrrDivergentFirstKmerBranchesFromPreviousBackbone) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTACGT", {});

    sharda::ReadPair pair;
    pair.read1.name = "orr_branch";
    pair.read1.seq = "TTTAC";
    pair.read1.cigar = {{sharda::CigarOp::M, 5}};
    pair.read1.ref_start = 1;
    pair.read1.ref_end = 6;
    pair.read1.flag = 0x43;

    pair.read2.name = "mate_unused";
    pair.read2.seq = "AC";
    pair.read2.cigar = {{sharda::CigarOp::M, 2}};
    pair.read2.ref_start = 0;
    pair.read2.ref_end = 2;
    pair.read2.flag = 0x83;

    sharda::add_read_pair(pair, graph, {});

    uint64_t first_read_node = graph.find_read_node("TTT");
    ASSERT_NE(first_read_node, UINT64_MAX);
    EXPECT_TRUE(graph.node(first_read_node).ref_positions.count(1) > 0);

    uint64_t predecessor = graph.backbone_node_at(0);
    ASSERT_NE(predecessor, UINT64_MAX);

    bool found_branch = false;
    for (uint64_t edge_index : graph.out_edges(predecessor)) {
        if (graph.edges()[edge_index].to == first_read_node) {
            found_branch = true;
            break;
        }
    }
    EXPECT_TRUE(found_branch);
}

TEST(ReadAdder, OrrReadNodesAccumulateImpliedCoordinates) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTACGT", {});

    sharda::ReadPair pair1;
    pair1.read1.name = "orr_coords_1";
    pair1.read1.seq = "TTTAC";
    pair1.read1.cigar = {{sharda::CigarOp::M, 5}};
    pair1.read1.ref_start = 1;
    pair1.read1.ref_end = 6;
    pair1.read1.flag = 0x43;
    pair1.read2.name = "mate_unused_1";
    pair1.read2.seq = "AC";
    pair1.read2.cigar = {{sharda::CigarOp::M, 2}};
    pair1.read2.ref_start = 0;
    pair1.read2.ref_end = 2;
    pair1.read2.flag = 0x83;

    sharda::ReadPair pair2;
    pair2.read1.name = "orr_coords_2";
    pair2.read1.seq = "TTTAC";
    pair2.read1.cigar = {{sharda::CigarOp::M, 5}};
    pair2.read1.ref_start = 4;
    pair2.read1.ref_end = 9;
    pair2.read1.flag = 0x43;
    pair2.read2.name = "mate_unused_2";
    pair2.read2.seq = "AC";
    pair2.read2.cigar = {{sharda::CigarOp::M, 2}};
    pair2.read2.ref_start = 0;
    pair2.read2.ref_end = 2;
    pair2.read2.flag = 0x83;

    sharda::add_read_pair(pair1, graph, {});
    sharda::add_read_pair(pair2, graph, {});

    uint64_t first_read_node = graph.find_read_node("TTT");
    ASSERT_NE(first_read_node, UINT64_MAX);
    EXPECT_TRUE(graph.node(first_read_node).ref_positions.count(1) > 0);
    EXPECT_TRUE(graph.node(first_read_node).ref_positions.count(4) > 0);
}

TEST_F(TempFileTest, DebugArtifactsWriteJsonAndViewer) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTAC", {});

    sharda::DebugArtifactsConfig config;
    config.enabled = true;
    config.output_dir = tmp_path("debug_artifacts");

    sharda::write_dbg_debug_artifacts(config, "raw", graph);

    EXPECT_TRUE(fs::exists(fs::path(config.output_dir) / "raw.json"));
    EXPECT_TRUE(fs::exists(fs::path(config.output_dir) / "raw.gfa"));
    EXPECT_TRUE(fs::exists(fs::path(config.output_dir) / "manifest.json"));
    EXPECT_TRUE(fs::exists(fs::path(config.output_dir) / "viewer.html"));

    std::ifstream in(fs::path(config.output_dir) / "raw.json");
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("\"graph_kind\": \"dbg\""), std::string::npos);
    EXPECT_NE(content.find("\"sequence\": \"ACG\""), std::string::npos);
}

TEST_F(TempFileTest, UnitigJsonIncludesNodeIds) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTAC", {});

    sharda::UnitigGraph ug;
    ASSERT_TRUE(ug.build(graph));

    std::string path = tmp_path("unitig.json");
    sharda::write_unitig_json(path, ug);

    std::ifstream in(path);
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("\"graph_kind\": \"unitig\""), std::string::npos);
    EXPECT_NE(content.find("\"node_ids\""), std::string::npos);
    EXPECT_NE(content.find("\"ref_pos\": 0"), std::string::npos);
    EXPECT_NE(content.find("\"ref_positions\": [0, 1, 2, 3]"), std::string::npos);
}

TEST_F(TempFileTest, UnitigGfaIncludesCoordinateTags) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTAC", {});

    sharda::UnitigGraph ug;
    ASSERT_TRUE(ug.build(graph));

    std::string path = tmp_path("unitig.gfa");
    sharda::write_unitig_gfa(path, ug);

    std::ifstream in(path);
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("\tRP:i:0\tRPS:Z:0,1,2,3"), std::string::npos);
}

namespace {

sharda::ReadPair make_trace_pair() {
    sharda::ReadPair pair;
    pair.read1.name = "trace_me";
    pair.read1.seq = "ACGTAC";
    pair.read1.cigar = {{sharda::CigarOp::M, 6}};
    pair.read1.ref_start = 0;
    pair.read1.ref_end = 6;
    pair.read1.flag = 0x43;
    pair.read2.name = "trace_me";
    pair.read2.seq = "AC";
    pair.read2.cigar = {{sharda::CigarOp::M, 2}};
    pair.read2.ref_start = 1;
    pair.read2.ref_end = 3;
    pair.read2.flag = 0x83;
    return pair;
}

std::vector<sharda::ReadTraceRecord> build_trace_records(sharda::DBG& graph) {
    auto pair = make_trace_pair();
    std::vector<std::string> traced_reads = {"trace_me"};
    std::vector<sharda::ReadTraceRecord> traces;
    sharda::ReadTraceSink trace_sink{&traced_reads, &traces};
    sharda::add_read_pair(pair, graph, {}, trace_sink);
    return traces;
}

} // namespace

TEST(ReadTrace, CapturesRawNodesBeforeUnitigBuild) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTAC", {});

    auto traces = build_trace_records(graph);

    ASSERT_EQ(traces.size(), 2u);
    ASSERT_FALSE(traces[0].raw_nodes.empty());
    EXPECT_TRUE(traces[1].raw_nodes.empty());
    EXPECT_EQ(traces[0].read_name, "trace_me");
}

TEST(ReadTrace, UnitigBuildSucceedsForTraceFixture) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTAC", {});

    auto traces = build_trace_records(graph);

    ASSERT_EQ(traces.size(), 2u);

    sharda::UnitigGraph ug;
    ASSERT_TRUE(ug.build(graph));
    EXPECT_NE(ug.node_to_unitig(traces[0].raw_nodes[0].node_id), UINT64_MAX);
}

TEST(ReadTrace, FinalizeReadTracesMapsUnitigs) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTAC", {});

    auto traces = build_trace_records(graph);

    sharda::UnitigGraph ug;
    ASSERT_TRUE(ug.build(graph));
    sharda::finalize_read_traces(traces, graph, ug);

    EXPECT_FALSE(traces[0].raw_nodes[0].removed_after_clean);
    EXPECT_NE(traces[0].raw_nodes[0].unitig_id, UINT64_MAX);
}

TEST_F(TempFileTest, ReadTraceArtifactsWriteJson) {
    sharda::DebugArtifactsConfig config;
    config.enabled = true;
    config.output_dir = tmp_path("debug_artifacts");
    config.traced_reads = {"trace_me"};

    std::vector<sharda::ReadTraceRecord> traces = {{
        "trace_me",
        "read1",
        "ORR",
        false,
        -1,
        {{0, "ACG", true, false, -1, -1, false, 0, "ACGT"}}
    }};

    sharda::write_read_trace_artifacts(config, traces);

    std::ifstream in(fs::path(config.output_dir) / "read_traces.json");
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("\"read_name\": \"trace_me\""), std::string::npos);
    EXPECT_NE(content.find("\"unitig_id\": 0"), std::string::npos);
}

TEST(LocusTrace, CapturesReferenceIntervalAndUnitigMapping) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTAC", {});

    auto traces = sharda::collect_locus_traces({{1, 3}}, "ACGTAC", 0, graph);
    ASSERT_EQ(traces.size(), 1u);
    EXPECT_EQ(traces[0].reference_sequence, "CGT");
    ASSERT_FALSE(traces[0].raw_nodes.empty());

    sharda::UnitigGraph ug;
    ASSERT_TRUE(ug.build(graph));
    sharda::finalize_locus_traces(traces, graph, ug);

    EXPECT_NE(traces[0].raw_nodes[0].unitig_id, UINT64_MAX);
}

TEST_F(TempFileTest, LocusTraceArtifactsWriteJson) {
    sharda::DebugArtifactsConfig config;
    config.enabled = true;
    config.output_dir = tmp_path("debug_artifacts");
    config.traced_loci = {{10, 5}};

    std::vector<sharda::LocusTraceRecord> traces = {{
        10,
        5,
        10,
        "ACGTA",
        {{0, "ACG", true, 10, -1, false, 0, "ACGT"}}
    }};

    sharda::write_locus_trace_artifacts(config, traces);

    std::ifstream in(fs::path(config.output_dir) / "locus_traces.json");
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("\"reference_sequence\": \"ACGTA\""), std::string::npos);
    EXPECT_NE(content.find("\"global_start\": 10"), std::string::npos);
}

TEST_F(TempFileTest, FlowPathArtifactsWriteJson) {
    sharda::DebugArtifactsConfig config;
    config.enabled = true;
    config.output_dir = tmp_path("debug_artifacts");

    std::vector<sharda::HaplotypePath> paths = {{
        {0, 4, 7},
        "ignored-sequence",
        12.5
    }};

    sharda::write_flow_path_artifacts(config, paths);

    std::ifstream in(fs::path(config.output_dir) / "flow_paths.json");
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("\"graph_kind\": \"flow_paths\""), std::string::npos);
    EXPECT_NE(content.find("\"flow\": 12.5"), std::string::npos);
    EXPECT_NE(content.find("\"unitig_ids\": [0, 4, 7]"), std::string::npos);
    EXPECT_EQ(content.find("\"sequence\""), std::string::npos);
    EXPECT_EQ(content.find("\"node_ids\""), std::string::npos);
}

TEST(DebugQuery, FindsReadTraceObject) {
        std::string text = R"JSON({
    "reads": [
        {"read_name": "read_a", "mate": "read1", "raw_nodes": []},
        {"read_name": "read_b", "mate": "read2", "raw_nodes": [{"node_id": 7}]}
    ]
})JSON";

        auto object = sharda::debug_query_find_read_trace_object(text, "read_b");
        ASSERT_TRUE(object.has_value());
        EXPECT_NE(object->find("\"read_name\": \"read_b\""), std::string::npos);
        EXPECT_NE(object->find("\"node_id\": 7"), std::string::npos);
}

TEST(DebugQuery, FindsLocusTraceObject) {
        std::string text = R"JSON({
    "loci": [
        {"local_start": 10, "length": 25, "reference_sequence": "AAAA", "raw_nodes": []},
        {"local_start": 100, "length": 50, "reference_sequence": "CCCC", "raw_nodes": [{"node_id": 4}]}
    ]
})JSON";

        auto object = sharda::debug_query_find_locus_trace_object(text, 100, 50);
        ASSERT_TRUE(object.has_value());
        EXPECT_NE(object->find("\"local_start\": 100"), std::string::npos);
        EXPECT_NE(object->find("\"node_id\": 4"), std::string::npos);
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

TEST(UnitigGraph, LinearCollapseWithReverseIdOrder) {
    sharda::DBG graph(3);
    auto sink = graph.add_backbone_node("GTA", 2, -1);
    auto mid = graph.add_backbone_node("CGT", 1, -1);
    auto source = graph.add_backbone_node("ACG", 0, -1);

    graph.add_edge(source, mid);
    graph.add_edge(mid, sink);

    sharda::UnitigGraph ug;
    bool ok = ug.build(graph);

    EXPECT_TRUE(ok);
    EXPECT_EQ(ug.unitig_count(), 1u);
    ASSERT_EQ(ug.unitigs().size(), 1u);
    EXPECT_EQ(ug.unitig(0).node_ids.size(), 3u);
}

TEST(UnitigGraph, DropsIsolatedSingletonKmerUnitigs) {
    sharda::DBG graph(3);
    sharda::build_backbone(graph, "ACGTAC", {}); // one connected unitig
    auto isolated = graph.add_read_node("TTT");

    sharda::UnitigGraph ug;
    bool ok = ug.build(graph);

    EXPECT_TRUE(ok);
    EXPECT_EQ(ug.unitig_count(), 1u);
    EXPECT_EQ(ug.node_to_unitig(isolated), UINT64_MAX);
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
