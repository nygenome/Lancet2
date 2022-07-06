#pragma once

#include "lancet2/core_enums.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class BaseCov {
 public:
  BaseCov() = default;

  u16 fwdRaw = 0;     // NOLINT
  u16 revRaw = 0;     // NOLINT
  u16 fwdBQPass = 0;  // NOLINT
  u16 revBQPass = 0;  // NOLINT

  [[nodiscard]] auto RawStrandCov(Strand s) const -> u16 { return s == Strand::FWD ? fwdRaw : revRaw; }
  [[nodiscard]] auto RawTotalCov() const -> u16 { return fwdRaw + revRaw; }

  [[nodiscard]] auto BQPassStrandCov(Strand s) const -> u16 { return s == Strand::FWD ? fwdBQPass : revBQPass; }

  [[nodiscard]] auto BQPassTotalCov() const -> u16 { return fwdBQPass + revBQPass; }
};
}  // namespace lancet2
