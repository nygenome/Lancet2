#pragma once

#include <string>

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
struct ContigInfo {
  std::string contigName;
  i64 contigLen;
};

inline auto CheckContigsMatch(absl::Span<const ContigInfo> hdr_ctgs, absl::Span<const ContigInfo> ref_ctgs) -> bool {
  absl::flat_hash_map<std::string, i64> refContigsMap;
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
}  // namespace lancet2
