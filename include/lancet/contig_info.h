#pragma once

#include <cstdint>
#include <string>

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"

namespace lancet {
struct ContigInfo {
  std::string contigName;
  std::int64_t contigLen;
};

inline auto CheckContigsMatch(absl::Span<const ContigInfo> hdr_ctgs, absl::Span<const ContigInfo> ref_ctgs) -> bool {
  absl::flat_hash_map<std::string, std::int64_t> refContigsMap;
  for (const auto& item : ref_ctgs) {
    refContigsMap[item.contigName] = item.contigLen;
  }

  for (const auto& ctg : hdr_ctgs) {
    if (!refContigsMap.contains(ctg.contigName) || refContigsMap.at(ctg.contigName) != ctg.contigLen) {
      return false;
    }
  }

  return true;
}
}  // namespace lancet
