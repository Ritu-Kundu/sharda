#pragma once

#include "graph/types.h"
#include <string>
#include <vector>

namespace sharda {

void write_vcf(const std::string& path,
               const std::vector<StructuralVariantCall>& calls,
               const std::string& source);

} // namespace sharda