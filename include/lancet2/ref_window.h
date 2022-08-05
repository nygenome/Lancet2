#pragma once

#include <string>
#include <string_view>
#include <utility>

#include "absl/strings/str_format.h"
#include "lancet2/genomic_region.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
// 0-based positions, includes start and excludes end
class RefWindow {
 public:
  RefWindow() = default;

  [[nodiscard]] auto GetWindowIndex() const -> usize { return winIdx; }
  [[nodiscard]] auto GetChromName() const -> std::string { return chromName; }
  [[nodiscard]] auto GetSequence() const -> std::string { return sequence; }
  [[nodiscard]] auto GetSeqView() const -> std::string_view { return sequence; }
  [[nodiscard]] auto GetStartPos0() const -> i64 { return startPos0; }
  [[nodiscard]] auto GetEndPos0() const -> i64 { return endPos0; }

  void SetWindowIndex(usize idx) { winIdx = idx; }
  void SetChromName(std::string chrom) { chromName = std::move(chrom); }
  void SetSequence(std::string seq) { sequence = std::move(seq); }
  void SetStartPos0(i64 start) { startPos0 = start; }
  void SetEndPos0(i64 end) { endPos0 = end; }

  [[nodiscard]] auto GetLength() const -> usize {
    if (startPos0 < 0 || endPos0 < 0) return 0;
    return static_cast<usize>(endPos0 - startPos0);
  }

  /// Convert region to samtools region string specification
  /// NOTE: samtools regions: 1-based positions, including start and end
  [[nodiscard]] auto ToRegionString() const -> std::string {
    if (startPos0 < 0 && endPos0 < 0) return chromName;
    if (startPos0 >= 0 && endPos0 < 0) return absl::StrFormat("%s:%d", chromName, startPos0 + 1);
    if (startPos0 < 0 && endPos0 >= 0) return absl::StrFormat("%s:-%d", chromName, endPos0);
    return absl::StrFormat("%s:%d-%d", chromName, startPos0 + 1, endPos0);
  }

  [[nodiscard]] auto ToSamtoolsRegion() const -> GenomicRegion {
    // resulting genomic region is 1-based, includes start and end
    return {chromName, static_cast<u32>(startPos0 + 1), static_cast<u32>(endPos0)};
  }

  [[nodiscard]] auto IsEmpty() const -> bool {
    return winIdx == 0 && chromName.empty() && sequence.empty() && startPos0 == -1 && endPos0 == -1;
  }

 private:
  usize winIdx = 0;
  std::string chromName;
  std::string sequence;

  i64 startPos0 = -1;
  i64 endPos0 = -1;
};
}  // namespace lancet2
