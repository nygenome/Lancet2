#include "lancet/base_hpcov.h"

#include <algorithm>
#include "lancet/assert_macro.h"

namespace lancet {
BaseHpCov::BaseHpCov(const BaseCov& cov) : BaseHpCov(cov, MakeDefaultHP(cov)) {}

BaseHpCov::BaseHpCov(const BaseCov& cov, const BaseHP& hp)
    : raw(std::make_pair(cov.fwdRaw, cov.revRaw), {hp[0].raw, hp[1].raw, hp[2].raw}),
      bqPass(std::make_pair(cov.fwdBQPass, cov.revBQPass), {hp[0].bqPass, hp[1].bqPass, hp[2].bqPass}) {}

auto BuildHPCovs(absl::Span<const BaseCov> covs) -> std::vector<BaseHpCov> {
  std::vector<BaseHpCov> result;
  result.reserve(covs.size());
  std::for_each(covs.cbegin(), covs.cend(), [&result](const auto& bcov) { result.emplace_back(bcov); });
  return result;
}

auto BuildHPCovs(absl::Span<const BaseCov> covs, absl::Span<const BaseHP> hps) -> std::vector<BaseHpCov> {
  LANCET_ASSERT(covs.size() == hps.size());  // NOLINT
  std::vector<BaseHpCov> result;
  result.reserve(covs.size());
  for (std::size_t idx = 0; idx < covs.size(); idx++) result.emplace_back(covs[idx], hps[idx]);
  return result;
}
}  // namespace lancet
