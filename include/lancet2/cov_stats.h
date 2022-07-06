#pragma once

#include <limits>

#include "lancet2/online_stats.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class CovStats {
 public:
  CovStats() = default;

  void Push(u16 val);

  [[nodiscard]] auto Mean() const -> float;
  [[nodiscard]] auto NonZeroMean() const -> float;

  [[nodiscard]] auto Minimum() const noexcept -> u16 { return stats.IsEmpty() ? 0 : allMin; }
  [[nodiscard]] auto NonZeroMinimum() const noexcept -> u16 { return non0Stats.IsEmpty() ? 0 : non0Min; }

 private:
  u16 allMin = std::numeric_limits<u16>::max();
  u16 non0Min = std::numeric_limits<u16>::max();

  OnlineStats stats;
  OnlineStats non0Stats;
};
}  // namespace lancet2
