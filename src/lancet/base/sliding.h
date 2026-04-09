#ifndef SRC_LANCET_BASE_SLIDING_H_
#define SRC_LANCET_BASE_SLIDING_H_

#include "lancet/base/assert.h"
#include "lancet/base/types.h"

#include "absl/container/fixed_array.h"
#include "absl/strings/string_view.h"

#include <string_view>

using SeqMers = absl::FixedArray<std::string_view>;
[[nodiscard]] inline auto SlidingView(std::string_view seq, usize const window) -> SeqMers {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (seq.length() < window)
    return absl::FixedArray<std::string_view>(0);

  auto const end_position = seq.length() - window;
  absl::FixedArray<std::string_view> result(end_position + 1);

  for (usize offset = 0; offset <= end_position; ++offset) {
    // NOLINTBEGIN(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
    result[offset] = absl::ClippedSubstr(seq, offset, window);
    LANCET_ASSERT(result[offset].length() == window)
    // NOLINTEND(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
  }

  return result;
}

#endif  // SRC_LANCET_BASE_SLIDING_H_
