#pragma once

#include <limits>

#include "lancet2/online_stats.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class CovStats {
 public:
  CovStats() = default;

  void Push(u16 val);

  [[nodiscard]] auto GetMean() const -> float;
  [[nodiscard]] auto GetNonZeroMean() const -> float;

  [[nodiscard]] auto GetMinimum() const noexcept -> u16 {
    return (stats.IsEmpty() || allMin == std::numeric_limits<u16>::max()) ? 0 : allMin;
  }

  [[nodiscard]] auto GetNonZeroMinimum() const noexcept -> u16 {
    return (non0Stats.IsEmpty() || non0Min == std::numeric_limits<u16>::max()) ? 0 : non0Min;
  }

 private:
  u16 allMin = std::numeric_limits<u16>::max();
  u16 non0Min = std::numeric_limits<u16>::max();

  OnlineStats stats;
  OnlineStats non0Stats;
};
}  // namespace lancet2
