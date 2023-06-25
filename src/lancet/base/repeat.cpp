#include "lancet/base/repeat.h"

#include <bit>
#include <memory>

#include "absl/container/flat_hash_set.h"
#include "lancet/base/assert.h"

// Based off of https://github.com/Daniel-Liu-c0deb0t/triple_accel/blob/master/src/hamming.rs
auto HammingDistWord64(std::string_view first, std::string_view second) -> usize {
  LANCET_ASSERT(first.length() == second.length())
  usize result = 0;

  const auto num_words = (first.length() >> 3);
  const auto rem_words = static_cast<unsigned long long>(first.length() & 7);

  // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
  const auto* aptr = reinterpret_cast<const unsigned long long*>(first.data());
  const auto* bptr = reinterpret_cast<const unsigned long long*>(second.data());
  // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)

  for (usize idx = 0; idx < num_words; ++idx) {
    auto val = (aptr[idx] ^ bptr[idx]);
    val |= val >> 4;
    val &= 0x0f0f0f0f0f0f0f0fULL;
    val |= val >> 2;
    val &= 0x3333333333333333ULL;
    val |= val >> 1;
    val &= 0x5555555555555555ULL;
    result += std::popcount(val);
  }

  if (rem_words > 0) {
    auto val = (aptr[num_words] ^ bptr[num_words]);
    val |= val >> 4;
    val &= 0x0f0f0f0f0f0f0f0fULL;
    val |= val >> 2;
    val &= 0x3333333333333333ULL;
    val |= val >> 1;
    val &= 0x5555555555555555ULL;
    // make sure to mask out bits outside the string lengths
    result += std::popcount((val & ((1ULL << ((rem_words) << 3ULL)) - 1ULL)));
  }

  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)

  return result;
}

auto HammingDistNaive(std::string_view first, std::string_view second) -> usize {
  LANCET_ASSERT(first.length() == second.length())
  usize result = 0;
  const auto length = first.length();
  for (usize idx = 0; idx < length; ++idx) {
    result += static_cast<usize>(first[idx] != second[idx]);
  }
  return result;
}

auto HasExactRepeat(absl::Span<const std::string_view> kmers) -> bool {
  const auto uniq_kmers = absl::flat_hash_set<std::string_view>(kmers.cbegin(), kmers.cend());
  return kmers.size() != uniq_kmers.size();
}

auto HasApproximateRepeat(absl::Span<const std::string_view> kmers, const i64 num_allowed_mismatches) -> bool {
  for (const auto& first_kmer : kmers) {
    for (const auto& second_kmer : kmers) {
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (std::addressof(first_kmer) == std::addressof(second_kmer)) continue;

      const auto dist = HammingDistWord64(first_kmer, second_kmer);
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (dist <= num_allowed_mismatches) return true;
    }
  }

  return false;
}
