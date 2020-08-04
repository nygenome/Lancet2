#pragma once

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

#include "absl/types/span.h"
#include "lancet/base_cov.h"
#include "lancet/base_hp.h"
#include "lancet/core_enums.h"

namespace lancet {
struct HpCov {
  HpCov(const std::pair<std::uint16_t, std::uint16_t>& cov, const std::array<std::uint16_t, 3>& hps)
      : fwdCov(cov.first), revCov(cov.second), HP0(hps[0]), HP1(hps[1]), HP2(hps[2]) {}
  HpCov() = delete;

  std::uint16_t fwdCov = 0;  // NOLINT
  std::uint16_t revCov = 0;  // NOLINT
  std::uint16_t HP0 = 0;     // NOLINT
  std::uint16_t HP1 = 0;     // NOLINT
  std::uint16_t HP2 = 0;     // NOLINT

  [[nodiscard]] auto TotalCov() const -> std::uint16_t { return fwdCov + revCov; }
};

struct BaseHpCov {
  explicit BaseHpCov(const BaseCov& cov);
  BaseHpCov(const BaseCov& cov, const BaseHP& hp);
  BaseHpCov() = delete;

  HpCov raw;     // NOLINT
  HpCov bqPass;  // NOLINT
};

[[nodiscard]] auto BuildHPCovs(absl::Span<const BaseCov> covs) -> std::vector<BaseHpCov>;
[[nodiscard]] auto BuildHPCovs(absl::Span<const BaseCov> covs, absl::Span<const BaseHP> hps) -> std::vector<BaseHpCov>;
}  // namespace lancet
