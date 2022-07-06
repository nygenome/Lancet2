#pragma once

#include <array>
#include <utility>
#include <vector>

#include "absl/types/span.h"
#include "lancet2/base_cov.h"
#include "lancet2/base_hp.h"
#include "lancet2/core_enums.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
struct HpCov {
  HpCov(const std::pair<u16, u16>& cov, const std::array<u16, 3>& hps)
      : fwdCov(cov.first), revCov(cov.second), HP0(hps[0]), HP1(hps[1]), HP2(hps[2]) {}
  HpCov() = delete;

  u16 fwdCov = 0;  // NOLINT
  u16 revCov = 0;  // NOLINT
  u16 HP0 = 0;     // NOLINT
  u16 HP1 = 0;     // NOLINT
  u16 HP2 = 0;     // NOLINT

  [[nodiscard]] auto GetTotalCov() const -> u16 { return fwdCov + revCov; }
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
}  // namespace lancet2
