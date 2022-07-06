#pragma once

#include <array>
#include <cstdint>

#include "lancet2/base_cov.h"

namespace lancet2 {
struct HPCount {
  std::uint16_t raw = 0;
  std::uint16_t bqPass = 0;
};

// HP0, HP1, HP2 -> HaplptypeCount
// HP0 -> Unassigned / Missing haplotype
// HP1 -> Haplotype 1
// HP2 -> Haplotype 2
using BaseHP = std::array<HPCount, 3>;

inline auto MakeDefaultHP(const BaseCov& cov) -> BaseHP {
  const auto baseRaw = cov.RawTotalCov();
  const auto baseBqp = cov.BQPassTotalCov();
  return BaseHP{HPCount{baseRaw, baseBqp}, HPCount{0, 0}, HPCount{0, 0}};
}
}  // namespace lancet2
