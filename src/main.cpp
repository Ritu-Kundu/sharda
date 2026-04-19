#include <iostream>
#include <string>
#include <thread>
#include <atomic>
#include <mutex>
#include <filesystem>
#include <spdlog/spdlog.h>

#include "util/log.h"
#include "io/fasta_reader.h"
#include "io/bed_reader.h"
#include "io/bam_reader.h"
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

namespace fs = std::filesystem;

namespace {

struct Args {
    std::string ref_fasta;
    std::string bam;
    std::string bed;          // TR BED
    std::string targets_bed;  // target regions BED (whole-genome mode)
    int         ploidy  = 2;
    int         k       = 121;
    int         threads = 1;
    int         padding = 1000;
    std::string out_prefix = "sharda_out";
    bool        debug  = false;
};

void usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " -r <ref.fa> -b <reads.bam> -p <ploidy>\n"
              << "       [-R <targets.bed>] [-j threads] [-f padding]\n"
              << "       [-t <repeats.bed>]\n"
              << "       [-k kmer_size] [-o out_prefix] [-d]\n"
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
              << "  -d  Debug logging\n";
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
        else if (arg == "-d") a.debug = true;
        else if (arg == "-h" || arg == "--help") { usage(argv[0]); std::exit(0); }
        else { std::cerr << "Unknown arg: " << arg << '\n'; usage(argv[0]); std::exit(1); }
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
    spdlog::info("Reading reference: {}", args.ref_fasta);
    auto [ref_name, ref_seq] = sharda::read_fasta(args.ref_fasta);
    spdlog::info("Reference: {} ({} bp)", ref_name, ref_seq.size());

    std::vector<sharda::TandemRepeat> trs;
    if (!args.bed.empty()) {
        spdlog::info("Reading tandem repeats: {}", args.bed);
        trs = sharda::read_bed(args.bed);
        spdlog::info("Tandem repeats: {}", trs.size());
    } else {
        spdlog::info("No tandem repeat BED provided; assembling without TR annotations");
    }

    // ── 2. Build backbone ───────────────────────────────────────────
    spdlog::info("Building backbone (k={})", args.k);
    sharda::DBG graph(args.k);
    sharda::build_backbone(graph, ref_seq, trs);

    // ── 3. Add reads ────────────────────────────────────────────────
    spdlog::info("Adding reads from: {}", args.bam);
    uint64_t read_pairs = 0;
    sharda::iterate_read_pairs(args.bam, [&](sharda::ReadPair&& pair) {
        sharda::add_read_pair(pair, graph, trs);
        read_pairs++;
    });
    spdlog::info("Added {} read pairs", read_pairs);
    spdlog::info("Graph after read addition: {} nodes, {} edges, {} hap_edges",
                 graph.node_count(), graph.edge_count(),
                 graph.haplotype_edges().size());

    std::string raw_gfa = args.out_prefix + ".raw.gfa";
    sharda::write_gfa(raw_gfa, graph);
    spdlog::info("Raw graph GFA: {}", raw_gfa);

    // ── 4. Clean graph ──────────────────────────────────────────────
    spdlog::info("Cleaning graph");
    int mean_read_len = 150;
    sharda::clean_graph(graph, mean_read_len);

    std::string clean_gfa = args.out_prefix + ".clean.gfa";
    sharda::write_gfa(clean_gfa, graph);
    spdlog::info("Cleaned graph GFA: {}", clean_gfa);

    // ── 5. Build unitig graph ───────────────────────────────────────
    spdlog::info("Building unitig graph");
    sharda::UnitigGraph ug;
    if (!ug.build(graph)) {
        spdlog::error("Unitig graph has cycles — aborting");
        return 1;
    }

    std::string unitig_gfa = args.out_prefix + ".unitig.gfa";
    sharda::write_unitig_gfa(unitig_gfa, ug);
    spdlog::info("Unitig graph GFA: {}", unitig_gfa);

    // ── 6. Flow decomposition ───────────────────────────────────────
    spdlog::info("Running flow decomposition (ploidy={})", args.ploidy);
    auto paths = sharda::flow_decomposition(ug, args.ploidy);

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
                params.debug        = args.debug;

                if (args.debug) {
                    std::string safe_name = region_name;
                    std::replace(safe_name.begin(), safe_name.end(), ':', '_');
                    params.debug_dir = (debug_base / safe_name).string();
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
