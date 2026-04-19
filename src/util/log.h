#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <memory>

namespace sharda {

/// Call once at startup to configure the global logger.
inline void init_logging(bool debug = false) {
    auto logger = spdlog::stdout_color_mt("sharda");
    logger->set_level(debug ? spdlog::level::debug : spdlog::level::info);
    spdlog::set_default_logger(logger);
    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
}

} // namespace sharda
