#include "lancet/fisher_exact.h"

#include <cmath>
#include <limits>

namespace lancet {
auto FisherExact::Test(int n_11, int n_12, int n_21, int n_22) -> FisherExact::Result {
  int i;         // NOLINT
  int j;         // NOLINT
  int max;       // NOLINT
  int min;       // NOLINT
  double p;      // NOLINT
  double q;      // NOLINT
  double left;   // NOLINT
  double right;  // NOLINT
  int nrt;       // NOLINT
  int nct;       // NOLINT
  int n;         // NOLINT
  FisherExact::Result result;
  HgaccT aux;

  nrt = n_11 + n_12;
  nct = n_11 + n_21;
  n = n_11 + n_12 + n_21 + n_22;  // calculate n1_, n_1 and n
  max = (nct < nrt) ? nct : nrt;  // max n11, for right tail
  min = nrt + nct - n;

  // min n11, for left tail
  if (min < 0) {
    min = 0;
  }

  // no need to do test
  if (min == max) {
    result.probability = 1.0;
    return result;
  }

  // the probability of the current table
  q = HypergeoAcc(n_11, nrt, nct, n, &aux);

  // left tail
  p = HypergeoAcc(min, 0, 0, 0, &aux);
  for (left = 0., i = min + 1; p < 0.99999999 * q; ++i) {  // NOLINT
    // loop until underflow
    left += p;
    p = HypergeoAcc(i, 0, 0, 0, &aux);
  }

  --i;
  if (p < 1.00000001 * q) {  // NOLINT
    left += p;
  } else {
    --i;
  }

  // right tail
  p = HypergeoAcc(max, 0, 0, 0, &aux);
  for (right = 0., j = max - 1; p < 0.99999999 * q; --j) {  // NOLINT
    // loop until underflow
    right += p;
    p = HypergeoAcc(j, 0, 0, 0, &aux);
  }

  ++j;
  if (p < 1.00000001 * q) {  // NOLINT
    right += p;
  } else {
    ++j;
  }

  // two-tail
  result.twoSidedPVal = left + right;
  if (result.twoSidedPVal > 1.0) result.twoSidedPVal = 1.0;

  // adjust left and right
  if (std::abs(i - n_11) < std::abs(j - n_11)) {
    right = 1. - left + q;
  } else {
    left = 1.0 - right + q;
  }

  result.leftPVal = left;
  result.rightPVal = right;
  result.probability = q;
  return result;
}

auto FisherExact::LBinom(int n, int k) -> double {
  if (k == 0 || n == k) {
    return 0;
  }
  return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
}

auto FisherExact::Hypergeo(int n_11, int nrt, int nct, int n) -> double {
  return std::exp(LBinom(nrt, n_11) + LBinom(n - nrt, nct - n_11) - LBinom(n, nct));
}

auto FisherExact::HypergeoAcc(int n_11, int nrt, int nct, int n, FisherExact::HgaccT *aux) -> double {
  if (static_cast<bool>(nrt) || static_cast<bool>(nct) || static_cast<bool>(n)) {
    aux->n11 = n_11;
    aux->n1_ = nrt;
    aux->n_1 = nct;
    aux->n = n;
  } else {  // then only n11 changed; the rest fixed
    if (static_cast<bool>(n_11 % 11) && static_cast<bool>(n_11 + aux->n - aux->n1_ - aux->n_1)) {  // NOLINT
      if (n_11 == aux->n11 + 1) {                                                                  // incremental
        aux->p *= static_cast<double>(aux->n1_ - aux->n11) / n_11 * (aux->n_1 - aux->n11) /
                  (n_11 + aux->n - aux->n1_ - aux->n_1);
        aux->n11 = n_11;
        return aux->p;
      }
      if (n_11 == aux->n11 - 1) {  // incremental
        aux->p *= static_cast<double>(aux->n11) / (aux->n1_ - n_11) * (aux->n11 + aux->n - aux->n1_ - aux->n_1) /
                  (aux->n_1 - n_11);
        aux->n11 = n_11;
        return aux->p;
      }
    }
    aux->n11 = n_11;
  }
  aux->p = Hypergeo(aux->n11, aux->n1_, aux->n_1, aux->n);
  return aux->p;
}

auto PhredScaledProbability(const FisherExact::Result &result) -> double {
  if (result.probability == 1.0) return 0.0;
  if (result.probability == 0.0) return -10.0 * std::log10(1 / std::numeric_limits<double>::max());  // NOLINT
  return -10.0 * std::log10(result.probability);                                                     // NOLINT
}
}  // namespace lancet
