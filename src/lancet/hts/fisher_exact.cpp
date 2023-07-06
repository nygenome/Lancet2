#include "lancet/hts/fisher_exact.h"

#include "htslib/kfunc.h"

namespace lancet::hts {

auto FisherExact::Test(const ContingencyTable &table) -> Result {
  const auto [case_counts, ctrl_counts] = table;
  const auto [n_11, n_12] = case_counts;
  const auto [n_21, n_22] = ctrl_counts;

  Result result;
  result.mDataProb = kt_fisher_exact(n_11, n_12, n_21, n_22, &result.mLessProb, &result.mMoreProb, &result.mDiffProb);
  return result;
}

}  // namespace lancet::hts
