#pragma once

#include "graph/types.h"
#include <string>
#include <vector>
#include <unordered_map>

namespace sharda {

class DBG; // forward declaration
class UnitigGraph;

/// Write a DBG to GFA1 format for debugging.
void write_gfa(const std::string& path, const DBG& graph);

/// Write a DBG to a JSON snapshot for structured debugging.
void write_dbg_json(const std::string& path, const DBG& graph);

/// Write a UnitigGraph to GFA1 format for debugging.
void write_unitig_gfa(const std::string& path, const UnitigGraph& graph);

/// Write an SV-oriented UnitigGraph GFA with support-class annotations.
void write_sv_unitig_gfa(const std::string& path, const UnitigGraph& graph);

/// Write a UnitigGraph to a JSON snapshot for structured debugging.
void write_unitig_json(const std::string& path, const UnitigGraph& graph);

/// Write an SV-oriented UnitigGraph JSON snapshot with support-class annotations.
void write_sv_unitig_json(const std::string& path, const UnitigGraph& graph);

} // namespace sharda
