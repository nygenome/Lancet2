#pragma once

#include <string>
#include <string_view>

#include "absl/container/fixed_array.h"
#include "lancet2/kmer.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
[[nodiscard]] auto KMovingSubstrs(std::string_view sv, usize k) -> absl::FixedArray<std::string>;

[[nodiscard]] auto CanonicalKmers(std::string_view sv, usize k) -> absl::FixedArray<Kmer>;

[[nodiscard]] auto CanonicalKmerHashes(std::string_view sv, usize k) -> absl::FixedArray<usize>;
}  // namespace lancet2
