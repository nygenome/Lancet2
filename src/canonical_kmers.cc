#include "lancet2/canonical_kmers.h"

#include "absl/strings/string_view.h"
#include "lancet2/assert_macro.h"

namespace lancet2 {
auto KMovingSubstrs(std::string_view sv, usize k) -> absl::FixedArray<std::string> {
  const auto endPos = sv.length() - k;
  absl::FixedArray<std::string> result(endPos + 1);

  for (usize offset = 0; offset <= endPos; offset++) {
    const auto subSeq = absl::ClippedSubstr(sv, offset, k);
    LANCET_ASSERT(subSeq.length() == k);  // NOLINT
    result[offset] = std::string(subSeq);
  }

  return result;
}

auto CanonicalKmers(std::string_view sv, usize k) -> absl::FixedArray<Kmer> {
  const auto endPos = sv.length() - k;
  absl::FixedArray<Kmer> result(endPos + 1);

  for (usize offset = 0; offset <= endPos; offset++) {
    const auto subSeq = absl::ClippedSubstr(sv, offset, k);
    LANCET_ASSERT(subSeq.length() == k);  // NOLINT
    result[offset] = Kmer(subSeq);
  }

  return result;
}

auto CanonicalKmerHashes(std::string_view sv, usize k) -> absl::FixedArray<usize> {
  const auto endPos = sv.length() - k;
  absl::FixedArray<usize> result(endPos + 1);

  for (usize offset = 0; offset <= endPos; offset++) {
    const auto subSeq = absl::ClippedSubstr(sv, offset, k);
    LANCET_ASSERT(subSeq.length() == k);  // NOLINT
    result[offset] = Kmer(subSeq).ID();
  }

  return result;
}
}  // namespace lancet2
