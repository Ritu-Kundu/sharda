#include <iostream>
#include <string>
#include <thread>
#include <atomic>
#include <mutex>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <vector>
#include <spdlog/spdlog.h>

#include "util/log.h"
#include "io/fasta_reader.h"
#include "io/bed_reader.h"
#include "io/bam_reader.h"
#include "io/debug_artifacts.h"
#include "io/fasta_writer.h"
#include "io/gfa_writer.h"
#include "graph/types.h"
#include "graph/dbg.h"
#include "graph/backbone.h"
#include "graph/unitig_graph.h"
#include "assembly/read_adder.h"
#include "assembly/graph_cleaner.h"
#include "assembly/flow_decomp.h"
#include "assembly/region_assembler.h"
#include "util/debug_config.h"
#include "util/debug_query.h"

namespace fs = std::filesystem;

namespace {

class StageTimer {
public:
    explicit StageTimer(std::string stage)
        : stage_(std::move(stage)), start_(std::chrono::steady_clock::now()) {
        spdlog::info("Stage start: {}", stage_);
    }

    void checkpoint(const char* message) const {
        spdlog::info("Stage progress: {} - {} ({} ms)",
                     stage_, message, elapsed_ms());
    }

    void finish() const {
        spdlog::info("Stage done: {} ({} ms)", stage_, elapsed_ms());
    }

private:
    long long elapsed_ms() const {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::steady_clock::now() - start_)
            .count();
    }

    std::string stage_;
    std::chrono::steady_clock::time_point start_;
};

struct Args {
    std::string ref_fasta;
    std::string bam;
    std::string bed;          // TR BED
    std::string targets_bed;  // target regions BED (whole-genome mode)
    std::string debug_dir;
    std::string debug_node;
    std::string debug_stage;
    std::string debug_read;
    std::string debug_locus;
    std::vector<std::string> trace_reads;
    std::vector<sharda::LocusTraceRequest> trace_loci;
    int         ploidy  = 2;
    int         k       = 121;
    int         threads = 1;
    int         padding = 1000;
    std::string out_prefix = "sharda_out";
    bool        stop_after_unitig_graph = false;
    bool        debug  = false;

    bool has_debug_node_lookup() const {
        return !debug_node.empty();
    }

    bool has_debug_read_lookup() const {
        return !debug_read.empty();
    }

    bool has_debug_locus_lookup() const {
        return !debug_locus.empty();
    }
};

bool parse_trace_locus_arg(const std::string& value, sharda::LocusTraceRequest& request) {
    size_t separator = value.find(':');
    if (separator == std::string::npos) {
        return false;
    }

    try {
        request.start = std::stoi(value.substr(0, separator));
        request.length = std::stoi(value.substr(separator + 1));
    } catch (const std::exception&) {
        return false;
    }

    return request.length >= 0;
}

sharda::DebugArtifactsConfig make_debug_artifacts_config(const Args& args) {
    sharda::DebugArtifactsConfig config;
    if (!args.debug) {
        return config;
    }

    config.enabled = true;
    config.output_dir = args.out_prefix + "_debug";
    config.traced_reads = args.trace_reads;
    config.traced_loci = args.trace_loci;
    return config;
}

void usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " -r <ref.fa> -b <reads.bam> -p <ploidy>\n"
              << "       [-R <targets.bed>] [-j threads] [-f padding]\n"
              << "       [-t <repeats.bed>]\n"
              << "       [-k kmer_size] [-o out_prefix] [--unitig-only] [-d] [--trace-read <name>] [--trace-locus <start:length>]\n"
              << "       [--debug-dir <dir> --debug-node <name> [--debug-stage <stage>]]\n"
              << "       [--debug-dir <dir> --debug-read <name>] [--debug-dir <dir> --debug-locus <start:length>]\n"
              << "\n"
              << "  -r  Reference FASTA (indexed .fai required for -R mode)\n"
              << "  -b  BAM file (name-sorted for single-region mode;\n"
              << "      coordinate-sorted + indexed for -R mode)\n"
              << "  -t  Tandem repeat BED (optional)\n"
              << "  -p  Ploidy\n"
              << "  -R  Target regions BED (enables whole-genome parallel mode)\n"
              << "  -j  Number of threads (default: 1, used with -R)\n"
              << "  -f  Flanking padding in bp (default: 1000, used with -R)\n"
              << "  -k  Kmer size (default: 121)\n"
              << "  -o  Output prefix (default: sharda_out)\n"
              << "  --unitig-only  Stop after unitig graph construction; skip ILP and haplotype FASTA output\n"
              << "  -d  Debug logging\n"
              << "  --trace-read   Trace a named read through raw, clean, and unitig stages\n"
              << "  --trace-locus  Trace a local reference interval through raw, clean, and unitig stages\n"
              << "  --debug-dir    Inspect an existing debug artifact directory\n"
              << "  --debug-node   Look up a node/segment name in debug GFA output\n"
              << "  --debug-read   Look up a traced read in read_traces.json\n"
              << "  --debug-locus  Look up a traced locus in locus_traces.json using start:length\n"
              << "  --debug-stage  Restrict lookup to raw, clean, or unitig\n";
}

std::vector<std::string> split_tab_fields(const std::string& line) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (start <= line.size()) {
        size_t end = line.find('\t', start);
        if (end == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, end - start));
        start = end + 1;
    }
    return fields;
}

std::string debug_stage_gfa_path(const std::string& debug_dir,
                                 const std::string& stage_name) {
    return (fs::path(debug_dir) / (stage_name + ".gfa")).string();
}

bool find_node_in_gfa(const std::string& gfa_path,
                      const std::string& node_name,
                      std::string& sequence,
                      std::vector<std::string>& tags) {
    std::ifstream in(gfa_path);
    if (!in) {
        return false;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] != 'S') {
            continue;
        }

        auto fields = split_tab_fields(line);
        if (fields.size() < 3 || fields[0] != "S") {
            continue;
        }
        if (fields[1] != node_name) {
            continue;
        }

        sequence = fields[2];
        tags.assign(fields.begin() + 3, fields.end());
        return true;
    }

    return false;
}

int run_debug_node_lookup(const Args& args) {
    if (args.debug_dir.empty() || args.debug_node.empty()) {
        std::cerr << "Error: --debug-dir and --debug-node must be provided together\n";
        return 1;
    }

    std::vector<std::string> stages;
    if (!args.debug_stage.empty()) {
        stages.push_back(args.debug_stage);
    } else {
        stages = {"raw", "clean", "unitig"};
    }

    for (const auto& stage_name : stages) {
        std::string gfa_path = debug_stage_gfa_path(args.debug_dir, stage_name);
        std::string sequence;
        std::vector<std::string> tags;
        if (!find_node_in_gfa(gfa_path, args.debug_node, sequence, tags)) {
            continue;
        }

        std::cout << "stage\t" << stage_name << '\n';
        std::cout << "node\t" << args.debug_node << '\n';
        std::cout << "sequence\t" << sequence << '\n';
        for (const auto& tag : tags) {
            std::cout << "tag\t" << tag << '\n';
        }
        return 0;
    }

    std::cerr << "Error: node " << args.debug_node << " not found in "
              << args.debug_dir;
    if (!args.debug_stage.empty()) {
        std::cerr << " for stage " << args.debug_stage;
    }
    std::cerr << '\n';
    return 1;
}

int run_debug_read_lookup(const Args& args) {
    if (args.debug_dir.empty() || args.debug_read.empty()) {
        std::cerr << "Error: --debug-dir and --debug-read must be provided together\n";
        return 1;
    }

    std::string artifact_path = (fs::path(args.debug_dir) / "read_traces.json").string();
    if (!fs::exists(artifact_path)) {
        std::cerr << "Error: read trace artifact not found: " << artifact_path << '\n'
                  << "Hint: rerun assembly with -d --trace-read READ_NAME to create read_traces.json\n";
        return 1;
    }
    auto object = sharda::debug_query_find_read_trace_object(
        sharda::debug_query_read_text_file(artifact_path), args.debug_read);
    if (!object) {
        std::cerr << "Error: read " << args.debug_read << " not found in " << artifact_path << '\n';
        return 1;
    }

    std::cout << *object << '\n';
    return 0;
}

int run_debug_locus_lookup(const Args& args) {
    if (args.debug_dir.empty() || args.debug_locus.empty()) {
        std::cerr << "Error: --debug-dir and --debug-locus must be provided together\n";
        return 1;
    }

    sharda::LocusTraceRequest request;
    if (!parse_trace_locus_arg(args.debug_locus, request)) {
        std::cerr << "Error: --debug-locus must be start:length\n";
        return 1;
    }

    std::string artifact_path = (fs::path(args.debug_dir) / "locus_traces.json").string();
    if (!fs::exists(artifact_path)) {
        std::cerr << "Error: locus trace artifact not found: " << artifact_path << '\n'
                  << "Hint: rerun assembly with -d --trace-locus START:LENGTH to create locus_traces.json\n";
        return 1;
    }
    auto object = sharda::debug_query_find_locus_trace_object(
        sharda::debug_query_read_text_file(artifact_path), request.start, request.length);
    if (!object) {
        std::cerr << "Error: locus " << args.debug_locus << " not found in " << artifact_path << '\n';
        return 1;
    }

    std::cout << *object << '\n';
    return 0;
}

Args parse_args(int argc, char* argv[]) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-r" && i + 1 < argc) a.ref_fasta = argv[++i];
        else if (arg == "-b" && i + 1 < argc) a.bam = argv[++i];
        else if (arg == "-t" && i + 1 < argc) a.bed = argv[++i];
        else if (arg == "-R" && i + 1 < argc) a.targets_bed = argv[++i];
        else if (arg == "-p" && i + 1 < argc) a.ploidy = std::stoi(argv[++i]);
        else if (arg == "-k" && i + 1 < argc) a.k = std::stoi(argv[++i]);
        else if (arg == "-j" && i + 1 < argc) a.threads = std::stoi(argv[++i]);
        else if (arg == "-f" && i + 1 < argc) a.padding = std::stoi(argv[++i]);
        else if (arg == "-o" && i + 1 < argc) a.out_prefix = argv[++i];
        else if (arg == "--unitig-only") a.stop_after_unitig_graph = true;
        else if (arg == "--trace-read" && i + 1 < argc) a.trace_reads.push_back(argv[++i]);
        else if (arg == "--trace-locus" && i + 1 < argc) {
            sharda::LocusTraceRequest request;
            if (!parse_trace_locus_arg(argv[++i], request)) {
                std::cerr << "Error: --trace-locus must be start:length\n";
                std::exit(1);
            }
            a.trace_loci.push_back(request);
        }
        else if (arg == "--debug-dir" && i + 1 < argc) a.debug_dir = argv[++i];
        else if (arg == "--debug-node" && i + 1 < argc) a.debug_node = argv[++i];
        else if (arg == "--debug-read" && i + 1 < argc) a.debug_read = argv[++i];
        else if (arg == "--debug-locus" && i + 1 < argc) a.debug_locus = argv[++i];
        else if (arg == "--debug-stage" && i + 1 < argc) a.debug_stage = argv[++i];
        else if (arg == "-d") a.debug = true;
        else if (arg == "-h" || arg == "--help") { usage(argv[0]); std::exit(0); }
        else { std::cerr << "Unknown arg: " << arg << '\n'; usage(argv[0]); std::exit(1); }
    }
    int debug_query_modes = 0;
    debug_query_modes += a.has_debug_node_lookup() ? 1 : 0;
    debug_query_modes += a.has_debug_read_lookup() ? 1 : 0;
    debug_query_modes += a.has_debug_locus_lookup() ? 1 : 0;
    if (debug_query_modes > 1) {
        std::cerr << "Error: choose only one of --debug-node, --debug-read, or --debug-locus\n";
        std::exit(1);
    }
    if (a.has_debug_node_lookup()) {
        if (a.debug_dir.empty() || a.debug_node.empty()) {
            std::cerr << "Error: --debug-dir and --debug-node must be used together\n";
            usage(argv[0]);
            std::exit(1);
        }
        if (!a.debug_stage.empty() && a.debug_stage != "raw"
            && a.debug_stage != "clean" && a.debug_stage != "unitig") {
            std::cerr << "Error: --debug-stage must be raw, clean, or unitig\n";
            std::exit(1);
        }
        return a;
    }
    if (a.has_debug_read_lookup()) {
        if (a.debug_dir.empty() || a.debug_read.empty()) {
            std::cerr << "Error: --debug-dir and --debug-read must be used together\n";
            usage(argv[0]);
            std::exit(1);
        }
        return a;
    }
    if (a.has_debug_locus_lookup()) {
        if (a.debug_dir.empty() || a.debug_locus.empty()) {
            std::cerr << "Error: --debug-dir and --debug-locus must be used together\n";
            usage(argv[0]);
            std::exit(1);
        }
        sharda::LocusTraceRequest request;
        if (!parse_trace_locus_arg(a.debug_locus, request)) {
            std::cerr << "Error: --debug-locus must be start:length\n";
            std::exit(1);
        }
        return a;
    }
    if ((!a.trace_reads.empty() || !a.trace_loci.empty()) && !a.debug) {
        std::cerr << "Error: --trace-read and --trace-locus require -d so artifacts can be written\n";
        std::exit(1);
    }
    if (a.ref_fasta.empty() || a.bam.empty()) {
        std::cerr << "Error: -r and -b are required\n";
        usage(argv[0]);
        std::exit(1);
    }
    if (a.threads < 1) a.threads = 1;
    if (a.padding < 0) a.padding = 0;
    return a;
}

/// Run the single-region pipeline (original behavior).
int run_single_region(const Args& args) {
    // ── 1. Read inputs ──────────────────────────────────────────────
    StageTimer total_timer("single-region pipeline");
    spdlog::info("Reading reference: {}", args.ref_fasta);
    StageTimer reference_timer("read reference fasta");
    auto [ref_name, ref_seq] = sharda::read_fasta(args.ref_fasta);
    reference_timer.finish();
    spdlog::info("Reference: {} ({} bp)", ref_name, ref_seq.size());

    std::vector<sharda::TandemRepeat> trs;
    if (!args.bed.empty()) {
        spdlog::info("Reading tandem repeats: {}", args.bed);
        StageTimer tr_timer("read tandem repeat bed");
        trs = sharda::read_bed(args.bed);
        tr_timer.finish();
        spdlog::info("Tandem repeats: {}", trs.size());
    } else {
        spdlog::info("No tandem repeat BED provided; assembling without TR annotations");
    }

    // ── 2. Build backbone ───────────────────────────────────────────
    spdlog::info("Building backbone (k={})", args.k);
    sharda::DBG graph(args.k);
    auto debug_artifacts = make_debug_artifacts_config(args);
    std::vector<sharda::ReadTraceRecord> read_traces;
    sharda::ReadTraceSink trace_sink{&debug_artifacts.traced_reads, &read_traces};
    std::vector<sharda::LocusTraceRecord> locus_traces;
    StageTimer backbone_timer("build backbone");
    sharda::build_backbone(graph, ref_seq, trs);
    backbone_timer.finish();

    // ── 3. Add reads ────────────────────────────────────────────────
    spdlog::info("Adding reads from: {}", args.bam);
    uint64_t read_pairs = 0;
    StageTimer read_timer("add reads");
    sharda::iterate_read_pairs(args.bam, [&](sharda::ReadPair&& pair) {
        sharda::add_read_pair(pair, graph, trs, trace_sink);
        read_pairs++;
    });
    read_timer.finish();
    spdlog::info("Added {} read pairs", read_pairs);
    spdlog::info("Graph after read addition: {} nodes, {} edges, {} hap_edges",
                 graph.node_count(), graph.edge_count(),
                 graph.haplotype_edges().size());

    if (debug_artifacts.should_trace_loci()) {
        locus_traces = sharda::collect_locus_traces(debug_artifacts.traced_loci,
                                                    ref_seq,
                                                    0,
                                                    graph);
    }

    if (debug_artifacts.should_write()) {
        sharda::write_dbg_debug_artifacts(debug_artifacts, "raw", graph);
        spdlog::info("Raw graph debug artifacts: {}", debug_artifacts.output_dir);
    } else {
        std::string raw_gfa = args.out_prefix + ".raw.gfa";
        sharda::write_gfa(raw_gfa, graph);
        spdlog::info("Raw graph GFA: {}", raw_gfa);
    }

    // ── 4. Clean graph ──────────────────────────────────────────────
    spdlog::info("Cleaning graph");
    int mean_read_len = 150;
    StageTimer clean_timer("clean graph");
    sharda::clean_graph(graph, mean_read_len);
    clean_timer.finish();

    if (debug_artifacts.should_write()) {
        sharda::write_dbg_debug_artifacts(debug_artifacts, "clean", graph);
        spdlog::info("Cleaned graph debug artifacts: {}", debug_artifacts.output_dir);
    } else {
        std::string clean_gfa = args.out_prefix + ".clean.gfa";
        sharda::write_gfa(clean_gfa, graph);
        spdlog::info("Cleaned graph GFA: {}", clean_gfa);
    }

    // ── 5. Build unitig graph ───────────────────────────────────────
    spdlog::info("Building unitig graph");
    sharda::UnitigGraph ug;
    StageTimer unitig_timer("build unitig graph");
    if (!ug.build(graph)) {
        spdlog::error("Unitig graph has cycles — aborting");
        return 1;
    }
    unitig_timer.finish();

    if (debug_artifacts.should_trace_reads()) {
        sharda::finalize_read_traces(read_traces, graph, ug);
        sharda::write_read_trace_artifacts(debug_artifacts, read_traces);
    }
    if (debug_artifacts.should_trace_loci()) {
        sharda::finalize_locus_traces(locus_traces, graph, ug);
        sharda::write_locus_trace_artifacts(debug_artifacts, locus_traces);
    }

    if (debug_artifacts.should_write()) {
        sharda::write_unitig_debug_artifacts(debug_artifacts, "unitig", ug);
        spdlog::info("Unitig graph debug artifacts: {}", debug_artifacts.output_dir);
    } else {
        std::string unitig_gfa = args.out_prefix + ".unitig.gfa";
        sharda::write_unitig_gfa(unitig_gfa, ug);
        spdlog::info("Unitig graph GFA: {}", unitig_gfa);
    }

    if (args.stop_after_unitig_graph) {
        spdlog::info("Stopping after unitig graph construction (--unitig-only)");
        total_timer.finish();
        return 0;
    }

    // ── 6. Flow decomposition ───────────────────────────────────────
    spdlog::info("Running flow decomposition (ploidy={})", args.ploidy);
    StageTimer flow_timer("flow decomposition");
    sharda::FlowBoundaryAnchors anchors;
    anchors.start_node_id = graph.backbone_node_at(0);
    if (static_cast<int>(ref_seq.size()) >= args.k) {
        anchors.end_node_id = graph.backbone_node_at(
            static_cast<int32_t>(ref_seq.size()) - args.k);
    }
    auto paths = sharda::flow_decomposition(ug, args.ploidy, anchors);
    flow_timer.finish();

    if (debug_artifacts.should_write()) {
        sharda::write_flow_path_artifacts(debug_artifacts, paths);
    }

    if (paths.empty()) {
        spdlog::warn("No haplotype paths found");
        return 1;
    }

    spdlog::info("Extracted {} haplotype paths", paths.size());

    // ── 7. Write output FASTA ───────────────────────────────────────
    std::vector<std::pair<std::string, std::string>> fasta_entries;
    for (size_t i = 0; i < paths.size(); ++i) {
        std::string name = ref_name + "_hap" + std::to_string(i + 1)
                         + "_flow" + std::to_string(static_cast<int>(paths[i].flow));
        fasta_entries.emplace_back(name, paths[i].sequence);
        spdlog::info("  Haplotype {}: {} bp, flow={:.1f}",
                     i + 1, paths[i].sequence.size(), paths[i].flow);
    }

    std::string out_fasta = args.out_prefix + ".haplotypes.fa";
    sharda::write_fasta(out_fasta, fasta_entries);
    spdlog::info("Output: {}", out_fasta);
    total_timer.finish();
    return 0;
}

/// Run the whole-genome parallel pipeline.
int run_whole_genome(const Args& args) {
    spdlog::info("Whole-genome mode: {} threads, {} bp padding",
                 args.threads, args.padding);

    // ── 1. Read target regions and TRs ──────────────────────────────
    spdlog::info("Reading target regions: {}", args.targets_bed);
    auto targets = sharda::read_target_regions(args.targets_bed);
    spdlog::info("Target regions: {}", targets.size());

    std::vector<sharda::TandemRepeat> all_trs;
    if (!args.bed.empty()) {
        spdlog::info("Reading tandem repeats: {}", args.bed);
        all_trs = sharda::read_bed(args.bed);
        spdlog::info("Tandem repeats: {}", all_trs.size());
    } else {
        spdlog::info("No tandem repeat BED provided; assembling without TR annotations");
    }

    // ── 2. Create temp directory for per-region BAMs ────────────────
    fs::path tmp_base = fs::temp_directory_path() / "sharda_tmp";
    fs::create_directories(tmp_base);

    // Debug output directory
    fs::path debug_base;
    if (args.debug) {
        debug_base = fs::path(args.out_prefix + "_debug");
        fs::create_directories(debug_base);
    }

    // ── 3. Assemble regions in parallel ─────────────────────────────
    std::vector<sharda::RegionResult> results(targets.size());
    std::atomic<size_t> next_region{0};
    std::atomic<size_t> completed{0};
    std::mutex log_mutex;

    auto worker = [&]() {
        while (true) {
            size_t idx = next_region.fetch_add(1, std::memory_order_relaxed);
            if (idx >= targets.size()) break;

            const auto& target = targets[idx];
            std::string region_name = target.chrom + ":"
                + std::to_string(target.start) + "-"
                + std::to_string(target.end);

            try {
                // Padded extraction window
                int32_t ext_start = std::max(static_cast<int32_t>(0),
                                             target.start - args.padding);
                int32_t ext_end   = target.end + args.padding;

                // Extract reference subsequence
                std::string ref_seq = sharda::read_fasta_region(
                    args.ref_fasta, target.chrom, ext_start, ext_end);

                // Create per-region name-sorted BAM
                std::string region_bam = (tmp_base
                    / (region_name + ".namesorted.bam")).string();
                // Replace colons/dashes in filename
                std::replace(region_bam.begin(), region_bam.end(), ':', '_');
                sharda::create_region_bam(
                    args.bam, target.chrom, ext_start, ext_end, region_bam);

                // Filter and adjust TRs to local coordinates
                auto local_trs = sharda::filter_trs_for_region(
                    all_trs, target, args.padding);

                // Set up region parameters
                sharda::RegionParams params;
                params.region_name  = region_name;
                params.ref_seq      = std::move(ref_seq);
                params.bam_path     = region_bam;
                params.coord_offset = ext_start;
                params.trs          = std::move(local_trs);
                params.ploidy       = args.ploidy;
                params.k            = args.k;
                params.stop_after_unitig_graph = args.stop_after_unitig_graph;
                params.debug        = args.debug;
                params.debug_artifacts = make_debug_artifacts_config(args);

                if (args.debug) {
                    std::string safe_name = region_name;
                    std::replace(safe_name.begin(), safe_name.end(), ':', '_');
                    params.debug_dir = (debug_base / safe_name).string();
                    params.debug_artifacts.output_dir = params.debug_dir;
                }

                // Assemble
                results[idx] = sharda::assemble_region(params);

                // Clean up temp BAM
                fs::remove(region_bam);

            } catch (const std::exception& e) {
                results[idx].region_name = region_name;
                results[idx].error = e.what();
                spdlog::error("[{}] Failed: {}", region_name, e.what());
            }

            size_t done = completed.fetch_add(1, std::memory_order_relaxed) + 1;
            spdlog::info("Progress: {}/{} regions completed", done, targets.size());
        }
    };

    // Launch worker threads
    int num_threads = std::min(args.threads, static_cast<int>(targets.size()));
    spdlog::info("Launching {} worker threads for {} regions",
                 num_threads, targets.size());

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back(worker);
    }
    for (auto& t : threads) {
        t.join();
    }

    // Clean up temp directory
    fs::remove_all(tmp_base);

    // ── 4. Collect and write results ────────────────────────────────
    std::vector<std::pair<std::string, std::string>> merged;
    int succeeded = 0, failed = 0;

    for (const auto& res : results) {
        if (res.success) {
            succeeded++;
            for (const auto& hap : res.haplotypes) {
                merged.push_back(hap);
            }

            // Per-region FASTA in debug folder
            if (args.debug && !debug_base.empty()) {
                std::string safe_name = res.region_name;
                std::replace(safe_name.begin(), safe_name.end(), ':', '_');
                std::string region_fasta =
                    (debug_base / safe_name / "haplotypes.fa").string();
                sharda::write_fasta(region_fasta, res.haplotypes);
            }
        } else {
            failed++;
            spdlog::warn("Region {} failed: {}", res.region_name, res.error);
        }
    }

    spdlog::info("Assembly complete: {}/{} regions succeeded, {} failed",
                 succeeded, targets.size(), failed);

    if (args.stop_after_unitig_graph) {
        if (succeeded == 0) {
            spdlog::error("No regions completed unitig graph construction");
            return 1;
        }
        spdlog::info("Stopped after unitig graph construction for all successful regions (--unitig-only)");
        return 0;
    }

    if (merged.empty()) {
        spdlog::error("No haplotypes assembled across any region");
        return 1;
    }

    std::string out_fasta = args.out_prefix + ".haplotypes.fa";
    sharda::write_fasta(out_fasta, merged);
    spdlog::info("Merged output: {} ({} haplotypes)", out_fasta, merged.size());

    return 0;
}

} // anonymous namespace

int main(int argc, char* argv[]) {
    auto args = parse_args(argc, argv);
    if (args.has_debug_node_lookup()) {
        return run_debug_node_lookup(args);
    }
    if (args.has_debug_read_lookup()) {
        return run_debug_read_lookup(args);
    }
    if (args.has_debug_locus_lookup()) {
        return run_debug_locus_lookup(args);
    }
    sharda::init_logging(args.debug);

    try {
        if (!args.targets_bed.empty()) {
            return run_whole_genome(args);
        } else {
            return run_single_region(args);
        }
    } catch (const std::exception& e) {
        spdlog::error("Fatal: {}", e.what());
        return 1;
    }

    spdlog::info("Done");
    return 0;
}
