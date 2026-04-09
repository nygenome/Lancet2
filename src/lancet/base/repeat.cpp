#include "lancet/base/repeat.h"

#include "lancet/base/assert.h"
#include "lancet/base/types.h"

#include "absl/container/flat_hash_set.h"
#include "absl/types/span.h"

#include <bit>
#include <memory>
#include <string_view>

// Based off of https://github.com/Daniel-Liu-c0deb0t/triple_accel/blob/master/src/hamming.rs
auto HammingDistWord64(std::string_view first, std::string_view second) -> usize {
  LANCET_ASSERT(first.length() == second.length())
  usize result = 0;

  auto const num_words = (first.length() >> 3);
  auto const rem_words = static_cast<unsigned long long>(first.length() & 7);

  // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* aptr = reinterpret_cast<unsigned long long const*>(first.data());
  auto const* bptr = reinterpret_cast<unsigned long long const*>(second.data());
  // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)

  for (usize idx = 0; idx < num_words; ++idx) {
    auto val = (aptr[idx] ^ bptr[idx]);
    val |= val >> 4;
    val &= 0x0f'0f'0f'0f'0f'0f'0f'0fULL;
    val |= val >> 2;
    val &= 0x33'33'33'33'33'33'33'33ULL;
    val |= val >> 1;
    val &= 0x55'55'55'55'55'55'55'55ULL;
    result += std::popcount(val);
  }

  if (rem_words > 0) {
    auto val = (aptr[num_words] ^ bptr[num_words]);
    val |= val >> 4;
    val &= 0x0f'0f'0f'0f'0f'0f'0f'0fULL;
    val |= val >> 2;
    val &= 0x33'33'33'33'33'33'33'33ULL;
    val |= val >> 1;
    val &= 0x55'55'55'55'55'55'55'55ULL;
    // make sure to mask out bits outside the string lengths
    result += std::popcount((val & ((1ULL << (rem_words << 3ULL)) - 1ULL)));
  }

  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)

  return result;
}

auto HammingDistNaive(std::string_view first, std::string_view second) -> usize {
  LANCET_ASSERT(first.length() == second.length())
  usize result = 0;
  auto const length = first.length();
  for (usize idx = 0; idx < length; ++idx) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
    result += static_cast<usize>(first[idx] != second[idx]);
  }
  return result;
}

auto HasExactRepeat(absl::Span<std::string_view const> kmers) -> bool {
  auto const uniq_kmers = absl::flat_hash_set<std::string_view>(kmers.cbegin(), kmers.cend());
  return kmers.size() != uniq_kmers.size();
}

auto HasApproximateRepeat(absl::Span<std::string_view const> kmers,
                          usize const num_allowed_mismatches) -> bool {
  for (auto const& first_kmer : kmers) {
    for (auto const& second_kmer : kmers) {
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (std::addressof(first_kmer) == std::addressof(second_kmer))
        continue;

      auto const dist = HammingDistWord64(first_kmer, second_kmer);
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (dist < num_allowed_mismatches)
        return true;
    }
  }

  return false;
}
