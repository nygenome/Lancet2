#include "lancet2/canonical_kmers.h"

#include "absl/strings/string_view.h"
#include "lancet2/assert_macro.h"

namespace lancet2 {
auto KMovingSubstrs(std::string_view sv, usize k) -> absl::FixedArray<std::string> {
  if (sv.length() < k) return absl::FixedArray<std::string>(0);

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
  if (sv.length() < k) return absl::FixedArray<Kmer>(0);

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
  if (sv.length() < k) return absl::FixedArray<usize>(0);

  const auto endPos = sv.length() - k;
  absl::FixedArray<usize> result(endPos + 1);

  for (usize offset = 0; offset <= endPos; offset++) {
    const auto subSeq = absl::ClippedSubstr(sv, offset, k);
    LANCET_ASSERT(subSeq.length() == k);  // NOLINT
    result[offset] = Kmer(subSeq).GetHash();
  }

  return result;
}
}  // namespace lancet2
