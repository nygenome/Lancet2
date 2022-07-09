#include "lancet2/fisher_exact.h"

#include <cmath>
#include <limits>

#include "htslib/kfunc.h"

namespace lancet2 {

auto FisherTest(int n_11, int n_12, int n_21, int n_22) -> FisherExactResult {
  FisherExactResult result;
  result.probability = kt_fisher_exact(n_11, n_12, n_21, n_22, &result.leftp, &result.rightp, &result.twosidep);
  return result;
}

auto PhredScaled(const FisherExactResult &result) -> double {
  if (result.probability == 1.0) return 0.0;
  if (result.probability == 0.0) return -10.0 * std::log10(1 / std::numeric_limits<double>::max());  // NOLINT
  return -10.0 * std::log10(result.probability);                                                     // NOLINT
}
}  // namespace lancet2
