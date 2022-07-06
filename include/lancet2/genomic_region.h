#pragma once

#include <string>
#include <utility>

#include "absl/strings/str_format.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
/// GenomicRegion is a 1-based fully closed interval that represents
///// a samtools region string, i.e. region includes start and end base
class GenomicRegion {
 public:
  explicit GenomicRegion(std::string chrom) : chromName(std::move(chrom)) {}
  GenomicRegion(std::string chrom, u32 start1, u32 end1)
      : chromName(std::move(chrom)), startPos1(static_cast<i64>(start1)), endPos1(static_cast<i64>(end1)) {}

  GenomicRegion() = delete;

  [[nodiscard]] auto GetChromName() const -> std::string { return chromName; }
  [[nodiscard]] auto GetStartPos1() const -> i64 { return startPos1; }
  [[nodiscard]] auto GetEndPos1() const -> i64 { return endPos1; }

  [[nodiscard]] auto GetLength() const -> usize {
    if (startPos1 <= 0 || endPos1 <= 0) return 0;
    return static_cast<usize>(endPos1 - startPos1 + 1);
  }

  /// Convert region to samtools region string specification
  /// NOTE: samtools regions: 1-based positions, including start and end
  [[nodiscard]] auto ToSamtoolsRegion() const -> std::string {
    if (startPos1 <= 0 && endPos1 <= 0) return chromName;
    if (startPos1 > 0 && endPos1 <= 0) return absl::StrFormat("%s:%d", chromName, startPos1);
    if (startPos1 <= 0 && endPos1 > 0) return absl::StrFormat("%s:-%d", chromName, endPos1);
    return absl::StrFormat("%s:%d-%d", chromName, startPos1, endPos1);
  }

 private:
  std::string chromName;
  i64 startPos1 = -1;
  i64 endPos1 = -1;
};
}  // namespace lancet2
