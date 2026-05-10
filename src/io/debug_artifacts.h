#pragma once

#include "assembly/flow_decomp.h"
#include "util/debug_config.h"
#include "util/locus_trace.h"
#include "util/read_trace.h"
#include <vector>
#include <string>

namespace sharda {

class DBG;
class UnitigGraph;

void ensure_debug_output_dir(const DebugArtifactsConfig& config);
void write_dbg_debug_artifacts(const DebugArtifactsConfig& config,
                               const std::string& stage_name,
                               const DBG& graph);
void write_unitig_debug_artifacts(const DebugArtifactsConfig& config,
                                  const std::string& stage_name,
                                  const UnitigGraph& graph);
std::vector<LocusTraceRecord> collect_locus_traces(
    const std::vector<LocusTraceRequest>& requests,
    const std::string& reference_sequence,
    int32_t coord_offset,
    const DBG& graph);
void finalize_locus_traces(std::vector<LocusTraceRecord>& traces,
                           const DBG& graph,
                           const UnitigGraph& unitig_graph);
void write_locus_trace_artifacts(const DebugArtifactsConfig& config,
                                 const std::vector<LocusTraceRecord>& traces);
void finalize_read_traces(std::vector<ReadTraceRecord>& traces,
                          const DBG& graph,
                          const UnitigGraph& unitig_graph);
void write_read_trace_artifacts(const DebugArtifactsConfig& config,
                                const std::vector<ReadTraceRecord>& traces);
void write_flow_path_artifacts(const DebugArtifactsConfig& config,
                               const std::vector<HaplotypePath>& paths);

} // namespace sharda