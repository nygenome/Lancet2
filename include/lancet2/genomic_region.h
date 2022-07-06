#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>

#include "absl/strings/str_format.h"

namespace lancet2 {
/// 1-based positions. includes start and end
class GenomicRegion {
 public:
  explicit GenomicRegion(std::string chrom) : chromName(std::move(chrom)) {}
  GenomicRegion(std::string chrom, std::int64_t start1, std::int64_t end1)
      : chromName(std::move(chrom)), startPos1(start1), endPos1(end1) {}
  GenomicRegion() = delete;

  [[nodiscard]] auto Chromosome() const -> std::string { return chromName; }
  [[nodiscard]] auto StartPosition1() const -> std::int64_t { return startPos1; }
  [[nodiscard]] auto EndPosition1() const -> std::int64_t { return endPos1; }

  [[nodiscard]] auto Length() const -> std::size_t {
    if (startPos1 <= 0 || endPos1 <= 0) return 0;
    return static_cast<std::size_t>(endPos1 - startPos1 + 1);
  }

  /// Convert region to samtools region string specification
  /// NOTE: samtools regions: 1-based positions, including start and end
  [[nodiscard]] auto ToRegionString() const -> std::string {
    if (startPos1 <= 0 && endPos1 <= 0) return chromName;
    if (startPos1 > 0 && endPos1 <= 0) return absl::StrFormat("%s:%d", chromName, startPos1);
    if (startPos1 <= 0 && endPos1 > 0) return absl::StrFormat("%s:-%d", chromName, endPos1);
    return absl::StrFormat("%s:%d-%d", chromName, startPos1, endPos1);
  }

 private:
  std::string chromName;
  std::int64_t startPos1 = -1;
  std::int64_t endPos1 = -1;
};
}  // namespace lancet2
