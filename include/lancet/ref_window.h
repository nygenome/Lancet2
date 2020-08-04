#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>

#include "absl/strings/str_format.h"
#include "lancet/genomic_region.h"

namespace lancet {
// 0-based positions, includes start and excludes end
class RefWindow {
 public:
  RefWindow() = default;

  [[nodiscard]] auto WindowIndex() const -> std::size_t { return winIdx; }
  [[nodiscard]] auto Chromosome() const -> std::string { return chromName; }
  [[nodiscard]] auto Sequence() const -> std::string { return sequence; }
  [[nodiscard]] auto SeqView() const -> std::string_view { return sequence; }
  [[nodiscard]] auto StartPosition0() const -> std::int64_t { return startPos0; }
  [[nodiscard]] auto EndPosition0() const -> std::int64_t { return endPos0; }

  void SetWindowIndex(std::size_t idx) { winIdx = idx; }
  void SetChromosome(std::string chrom) { chromName = std::move(chrom); }
  void SetSequence(std::string seq) { sequence = std::move(seq); }
  void SetStartPosition0(std::int64_t start) { startPos0 = start; }
  void SetEndPosition0(std::int64_t end) { endPos0 = end; }

  [[nodiscard]] auto Length() const -> std::size_t {
    if (startPos0 < 0 || endPos0 < 0) return 0;
    return static_cast<std::size_t>(endPos0 - startPos0);
  }

  /// Convert region to samtools region string specification
  /// NOTE: samtools regions: 1-based positions, including start and end
  [[nodiscard]] auto ToRegionString() const -> std::string {
    if (startPos0 < 0 && endPos0 < 0) return chromName;
    if (startPos0 >= 0 && endPos0 < 0) return absl::StrFormat("%s:%d", chromName, startPos0 + 1);
    if (startPos0 < 0 && endPos0 >= 0) return absl::StrFormat("%s:-%d", chromName, endPos0);
    return absl::StrFormat("%s:%d-%d", chromName, startPos0 + 1, endPos0);
  }

  [[nodiscard]] auto ToGenomicRegion() const -> GenomicRegion {
    // resulting genomic region is 1-based, includes start and end
    return GenomicRegion(chromName, startPos0 + 1, endPos0);
  }

 private:
  std::size_t winIdx = 0;
  std::string chromName;
  std::string sequence;

  std::int64_t startPos0 = -1;
  std::int64_t endPos0 = -1;
};
}  // namespace lancet
