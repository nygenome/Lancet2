#include "lancet2/cov_stats.h"

#include <algorithm>

namespace lancet2 {
void CovStats::Push(u16 val) {
  allMin = std::min(val, allMin);
  stats.Add(val);

  if (val > 0) {
    non0Min = std::min(val, non0Min);
    non0Stats.Add(val);
  }
}

auto CovStats::GetMean() const -> float { return static_cast<float>(stats.GetMean()); }
auto CovStats::GetNonZeroMean() const -> float { return static_cast<float>(non0Stats.GetMean()); }
}  // namespace lancet2
