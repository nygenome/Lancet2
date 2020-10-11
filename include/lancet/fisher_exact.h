#pragma once

namespace lancet {
struct FisherExactResult {
  double leftp = 0.0;
  double rightp = 0.0;
  double twosidep = 0.0;
  double probability = 0.0;
};

[[nodiscard]] auto FisherTest(int n_11, int n_12, int n_21, int n_22) -> FisherExactResult;
[[nodiscard]] auto PhredScaled(const FisherExactResult& result) -> double;

[[nodiscard]] inline auto PhredFisherScore(int n_11, int n_12, int n_21, int n_22) -> double {
  return PhredScaled(FisherTest(n_11, n_12, n_21, n_22));
}
}  // namespace lancet
