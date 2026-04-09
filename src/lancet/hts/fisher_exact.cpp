#include "lancet/hts/fisher_exact.h"

#include "htslib/kfunc.h"

namespace lancet::hts {

auto FisherExact::Test(ContingencyTable const& table) -> Result {
  auto const [case_counts, ctrl_counts] = table;
  auto const [n_11, n_12] = case_counts;
  auto const [n_21, n_22] = ctrl_counts;

  Result result;
  result.mDataProb = kt_fisher_exact(n_11, n_12, n_21, n_22, &result.mLessProb, &result.mMoreProb,
                                     &result.mDiffProb);
  return result;
}

}  // namespace lancet::hts
