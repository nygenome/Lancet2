#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "absl/types/span.h"
#include "lancet2/core_enums.h"

namespace lancet2 {
[[nodiscard]] auto CanMergeSeqs(std::string_view source, std::string_view buddy, BuddyPosition merge_dir,
                                bool reverse_buddy, std::size_t k) -> bool;

void MergeKmerSeqs(std::string* source, std::string_view buddy, BuddyPosition dir, bool reverse_buddy, std::size_t k);

template <typename T>
void MergeNodeInfo(std::vector<T>* source, absl::Span<const T> buddy, BuddyPosition dir, bool reverse_buddy,
                   std::size_t k) {
  assert(source != nullptr);

  if (source->empty()) {
    reverse_buddy ? source->insert(source->end(), buddy.crbegin(), buddy.crend())
                  : source->insert(source->end(), buddy.cbegin(), buddy.cend());
    return;
  }

  switch (dir) {
    case BuddyPosition::FRONT:
      reverse_buddy ? source->insert(source->end(), buddy.crbegin() + k - 1, buddy.crend())
                    : source->insert(source->end(), buddy.cbegin() + k - 1, buddy.cend());
      break;
    case BuddyPosition::BACK:
    default:
      reverse_buddy ? source->insert(source->begin(), buddy.crbegin(), buddy.crend() - k + 1)
                    : source->insert(source->begin(), buddy.cbegin(), buddy.cend() - k + 1);
      break;
  }
}
}  // namespace lancet2
