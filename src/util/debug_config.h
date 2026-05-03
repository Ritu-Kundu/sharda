#pragma once

#include <string>

namespace sharda {

struct DebugArtifactsConfig {
    bool        enabled        = false;
    std::string output_dir;
    bool        emit_gfa       = true;
    bool        emit_json      = true;
    bool        emit_html_view = true;

    bool should_write() const {
        return enabled && !output_dir.empty();
    }
};

} // namespace sharda