#pragma once

#include <cstddef>
#include <string>
#include <string_view>

#include "absl/container/fixed_array.h"
#include "lancet2/kmer.h"

namespace lancet2 {
[[nodiscard]] auto KMovingSubstrs(std::string_view sv, std::size_t k) -> absl::FixedArray<std::string>;

[[nodiscard]] auto CanonicalKmers(std::string_view sv, std::size_t k) -> absl::FixedArray<Kmer>;

[[nodiscard]] auto CanonicalKmerHashes(std::string_view sv, std::size_t k) -> absl::FixedArray<std::size_t>;
}  // namespace lancet2
