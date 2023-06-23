#ifndef SRC_LANCET_BASE_REPEAT_H_
#define SRC_LANCET_BASE_REPEAT_H_

#include <string_view>

#include "absl/types/span.h"
#include "lancet/base/types.h"

[[nodiscard]] auto HammingDistWord64(std::string_view first, std::string_view second) -> usize;
[[nodiscard]] auto HammingDistNaive(std::string_view first, std::string_view second) -> usize;

[[nodiscard]] auto HasExactRepeat(absl::Span<const std::string_view> kmers) -> bool;
[[nodiscard]] auto HasApproximateRepeat(absl::Span<const std::string_view> kmers, i64 num_allowed_mismatches) -> bool;

#endif  // SRC_LANCET_BASE_REPEAT_H_
