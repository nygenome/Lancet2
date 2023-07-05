#ifndef SRC_LANCET_BASE_SLIDING_H_
#define SRC_LANCET_BASE_SLIDING_H_

#include <string_view>

#include "absl/container/fixed_array.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
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

template <typename T>
using DataMers = absl::FixedArray<absl::Span<const T>>;
template <typename T>
[[nodiscard]] inline auto SlidingView(absl::Span<const T> data, const usize window) -> DataMers<T> {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (data.length() < window) return absl::FixedArray<absl::Span<const T>>(0);

  const auto end_position = data.length() - window;
  absl::FixedArray<absl::Span<const T>> result(end_position + 1);

  for (usize offset = 0; offset <= end_position; ++offset) {
    result[offset] = data.subspan(offset, window);
    LANCET_ASSERT(result[offset].length() == window)
  }

  return result;
}

#endif  // SRC_LANCET_BASE_SLIDING_H_
