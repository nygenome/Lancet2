#ifndef SRC_LANCET_HTS_FISHER_EXACT_H_
#define SRC_LANCET_HTS_FISHER_EXACT_H_

#include <array>

#include "lancet/base/types.h"

namespace lancet::hts {

/// Fisher's Exact Test is a statistical significance test used to determine
/// if there are nonrandom associations between two categorical variables,
/// organized in a 2x2 contingency table.
///
/// It calculates the probability distribution of data as extreme (or more)
/// under the null hypothesis that the variables are independent of each other.
///
/// To clarify, let's use an example of a 2x2 contingency table for a study
/// where we are comparing cases and controls:
///
///  |      |  variant   | non-variant |
///  |------|------------|-------------|
///  | Case |    n_11    |     n_12    |
///  | Ctrl |    n_21    |     n_22    |
///
/// `n_11` is the number of cases with the variant
/// `n_12` is the number of cases without the variant
/// `n_21` is the number of controls with the variant
/// `n_22` is the number of controls without the variant
class FisherExact {
 public:
  /// The result of Fisher's Exact Test in this case might include:
  ///
  /// *  `mLessProb` (or `left`): The probability of obtaining a result as extreme as, or more
  ///     extreme than, the observed data, assuming that the true ratio of odds is less than 1.
  ///     This tests the null hypothesis that the odds ratio is >= 1.
  ///     This is useful when we have a specific directional hypothesis,
  ///     e.g., we believe the variant is less likely to occur in the cases than controls.
  ///
  /// *  `mMoreProb` (or `right`): The probability of obtaining a result as extreme as, or more
  ///     extreme than, the observed data, assuming that the true ratio of odds is greater than 1.
  ///     This tests the null hypothesis that the odds ratio is <= 1.
  ///     This is useful when we have a specific directional hypothesis,
  ///     e.g., we believe the variant is more likely to occur in the cases than controls.
  ///
  /// *  `mDiffProb` (or `two.sided`): The probability of obtaining a result as extreme as, or more
  ///     extreme than, the observed data in either direction, regardless of the direction of the true odds ratio.
  ///     This tests the null hypothesis that the odds ratio is 1.
  ///     This is useful when we don't have a specific directional hypothesis and just want to know
  ///     if there is a difference in the occurrence of the variant between the cases and controls.
  ///
  /// *  `mDataProb`: The probability of the observed data given the null hypothesis.
  ///     This is not typically used to calculate a p-value, but could give you an idea of
  ///     how likely your observed data is assuming the null hypothesis is true.
  ///
  ///  When calculating phred-scaled values (log-scaled base 10), we generally use the two-tailed p-value
  ///  to convert to the phred scale because it gives us an overall measure of statistical significance
  ///  without assuming a particular direction of the effect. However, if we have a specific hypothesis
  ///  about the direction of the effect, we could use the left or right p-values as appropriate.
  struct MidpAdjustedResult {
    f64 mLessProb = 0.0;
    f64 mMoreProb = 0.0;
    f64 mDiffProb = 0.0;
    f64 mDataProb = 0.0;
  };

  using Row = std::array<int, 2>;
  using ContingencyTable = std::array<Row, 2>;
  [[nodiscard]] static auto Test(const ContingencyTable& table) -> MidpAdjustedResult;
};

static constexpr u8 MAX_PHRED_SCORE = 255;
[[nodiscard]] auto PhredToErrorProb(u32 phred_score) -> f64;
[[nodiscard]] auto ErrorProbToPhred(f64 prob) -> u8;
[[nodiscard]] auto ClampPhredScore(f64 score) -> u8;

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_FISHER_EXACT_H_
