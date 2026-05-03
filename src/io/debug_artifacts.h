#pragma once

#include "util/debug_config.h"
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

} // namespace sharda