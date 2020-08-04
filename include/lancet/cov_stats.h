#pragma once

#include <cstdint>
#include <limits>

#include "lancet/online_stats.h"

namespace lancet {
class CovStats {
 public:
  CovStats() = default;

  void Push(std::uint16_t val);

  [[nodiscard]] auto Mean() const -> float;
  [[nodiscard]] auto NonZeroMean() const -> float;

  [[nodiscard]] auto Minimum() const noexcept -> std::uint16_t { return stats.IsEmpty() ? 0 : allMin; }
  [[nodiscard]] auto NonZeroMinimum() const noexcept -> std::uint16_t { return non0Stats.IsEmpty() ? 0 : non0Min; }

 private:
  std::uint16_t allMin = std::numeric_limits<std::uint16_t>::max();
  std::uint16_t non0Min = std::numeric_limits<std::uint16_t>::max();

  OnlineStats stats;
  OnlineStats non0Stats;
};
}  // namespace lancet
