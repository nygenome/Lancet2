#include "lancet2/canonical_kmers.h"

#include "lancet2/assert_macro.h"

#include "absl/strings/string_view.h"

namespace lancet2 {
auto KMovingSubstrs(std::string_view sv, std::size_t k) -> absl::FixedArray<std::string> {
  const auto endPos = sv.length() - k;
  absl::FixedArray<std::string> result(endPos + 1);

  for (std::size_t offset = 0; offset <= endPos; offset++) {
    const auto subSeq = absl::ClippedSubstr(sv, offset, k);
    LANCET_ASSERT(subSeq.length() == k);  // NOLINT
    result[offset] = std::string(subSeq);
  }

  return result;
}

auto CanonicalKmers(std::string_view sv, std::size_t k) -> absl::FixedArray<Kmer> {
  const auto endPos = sv.length() - k;
  absl::FixedArray<Kmer> result(endPos + 1);

  for (std::size_t offset = 0; offset <= endPos; offset++) {
    const auto subSeq = absl::ClippedSubstr(sv, offset, k);
    LANCET_ASSERT(subSeq.length() == k);  // NOLINT
    result[offset] = Kmer(subSeq);
  }

  return result;
}

auto CanonicalKmerHashes(std::string_view sv, std::size_t k) -> absl::FixedArray<std::size_t> {
  const auto endPos = sv.length() - k;
  absl::FixedArray<std::size_t> result(endPos + 1);

  for (std::size_t offset = 0; offset <= endPos; offset++) {
    const auto subSeq = absl::ClippedSubstr(sv, offset, k);
    LANCET_ASSERT(subSeq.length() == k);  // NOLINT
    result[offset] = Kmer(subSeq).ID();
  }

  return result;
}
}  // namespace lancet2
