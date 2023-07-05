#ifndef SRC_LANCET_BASE_SLIDING_H_
#define SRC_LANCET_BASE_SLIDING_H_

#include <string_view>

#include "absl/container/fixed_array.h"
#include "absl/strings/string_view.h"
#include "lancet/base/assert.h"
#include "lancet/base/types.h"

using SeqMers = absl::FixedArray<std::string_view>;
[[nodiscard]] inline auto SlidingView(std::string_view seq, const usize window) -> SeqMers {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (seq.length() < window) return absl::FixedArray<std::string_view>(0);

  const auto end_position = seq.length() - window;
  absl::FixedArray<std::string_view> result(end_position + 1);

  for (usize offset = 0; offset <= end_position; ++offset) {
    result[offset] = absl::ClippedSubstr(seq, offset, window);
    LANCET_ASSERT(result[offset].length() == window)
  }

  return result;
}

#endif  // SRC_LANCET_BASE_SLIDING_H_
