#pragma once

#include <cstdint>

#include "lancet2/core_enums.h"

namespace lancet2 {
class BaseCov {
 public:
  BaseCov() = default;

  std::uint16_t fwdRaw = 0;     // NOLINT
  std::uint16_t revRaw = 0;     // NOLINT
  std::uint16_t fwdBQPass = 0;  // NOLINT
  std::uint16_t revBQPass = 0;  // NOLINT

  [[nodiscard]] auto RawStrandCov(Strand s) const -> std::uint16_t { return s == Strand::FWD ? fwdRaw : revRaw; }
  [[nodiscard]] auto RawTotalCov() const -> std::uint16_t { return fwdRaw + revRaw; }

  [[nodiscard]] auto BQPassStrandCov(Strand s) const -> std::uint16_t {
    return s == Strand::FWD ? fwdBQPass : revBQPass;
  }

  [[nodiscard]] auto BQPassTotalCov() const -> std::uint16_t { return fwdBQPass + revBQPass; }
};
}  // namespace lancet2
