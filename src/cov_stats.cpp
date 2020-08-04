#include "lancet/cov_stats.h"

#include <algorithm>

namespace lancet {
void CovStats::Push(std::uint16_t val) {
  allMin = std::min(val, allMin);
  stats.Add(val);

  if (val > 0) {
    non0Min = std::min(val, non0Min);
    non0Stats.Add(val);
  }
}

auto CovStats::Mean() const -> float { return static_cast<float>(stats.Mean()); }
auto CovStats::NonZeroMean() const -> float { return static_cast<float>(non0Stats.Mean()); }
}  // namespace lancet
