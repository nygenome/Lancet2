#ifndef SRC_LANCET_BASE_REPEAT_H_
#define SRC_LANCET_BASE_REPEAT_H_

#include "lancet/base/types.h"

#include "absl/types/span.h"

#include <string_view>

[[nodiscard]] auto HammingDistWord64(std::string_view first, std::string_view second) -> usize;
[[nodiscard]] auto HammingDistNaive(std::string_view first, std::string_view second) -> usize;

[[nodiscard]] auto HasExactRepeat(absl::Span<std::string_view const> kmers) -> bool;
[[nodiscard]] auto HasApproximateRepeat(absl::Span<std::string_view const> kmers,
                                        usize num_allowed_mismatches) -> bool;

#endif  // SRC_LANCET_BASE_REPEAT_H_
