#include "lancet/base/mann_whitney.h"

#include "lancet/base/types.h"

#include "absl/types/span.h"
#include "catch_amalgamated.hpp"

#include <vector>

#include <cmath>

namespace lancet::base::tests {

// ============================================================================
// MannWhitneyEffectSize: NaN Remediation Tests
//
// Validates the std::optional return path:
//   nullopt → one or both groups empty (untestable)
//   0.0     → test ran, all values identical (genuine zero)
//   value   → test ran, found bias
// ============================================================================

TEST_CASE("MannWhitneyEffectSize returns nullopt when REF group is empty",
          "[lancet][base][MannWhitneyEffectSize]") {
  std::vector<f64> const ref_vals = {};
  std::vector<f64> const alt_vals = {1.0, 2.0, 3.0};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE_FALSE(result.has_value());
}

TEST_CASE("MannWhitneyEffectSize returns nullopt when ALT group is empty",
          "[lancet][base][MannWhitneyEffectSize]") {
  std::vector<f64> const ref_vals = {1.0, 2.0, 3.0};
  std::vector<f64> const alt_vals = {};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE_FALSE(result.has_value());
}

TEST_CASE("MannWhitneyEffectSize returns nullopt when both groups are empty",
          "[lancet][base][MannWhitneyEffectSize]") {
  std::vector<f64> const ref_vals = {};
  std::vector<f64> const alt_vals = {};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE_FALSE(result.has_value());
}

TEST_CASE("MannWhitneyEffectSize returns 0.0 when all values are identical",
          "[lancet][base][MannWhitneyEffectSize]") {
  std::vector<f64> const ref_vals = {5.0, 5.0, 5.0};
  std::vector<f64> const alt_vals = {5.0, 5.0, 5.0};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE(result.has_value());
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  REQUIRE(result.value() == Catch::Approx(0.0).margin(1e-10));
}

TEST_CASE("MannWhitneyEffectSize returns nonzero for biased groups",
          "[lancet][base][MannWhitneyEffectSize]") {
  // ALT values are clearly higher than REF → positive effect size
  std::vector<f64> const ref_vals = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<f64> const alt_vals = {6.0, 7.0, 8.0, 9.0, 10.0};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE(result.has_value());
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  REQUIRE(result.value() > 0.0);
}

TEST_CASE("MannWhitneyEffectSize works with u8 type", "[lancet][base][MannWhitneyEffectSize]") {
  std::vector<u8> const ref_vals = {60, 60, 60};
  std::vector<u8> const alt_vals = {};
  auto const result =
      MannWhitneyEffectSize<u8>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE_FALSE(result.has_value());
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  Sample-size edge cases                                                  ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("MannWhitneyEffectSize handles a single-element ALT group",
          "[lancet][base][MannWhitneyEffectSize]") {
  // Smallest non-empty ALT group. The variance formula divides by N*(N-1),
  // so the smallest legal N is 2 (1 ref + 1 alt) — covered below; this case
  // exercises N = 4 (3 ref + 1 alt) where the alt rank sum has just one
  // contribution.
  std::vector<f64> const ref_vals = {1.0, 2.0, 3.0};
  std::vector<f64> const alt_vals = {10.0};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE(result.has_value());
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  CHECK(result.value() > 0.0);
}

TEST_CASE("MannWhitneyEffectSize handles a single-element REF group",
          "[lancet][base][MannWhitneyEffectSize]") {
  std::vector<f64> const ref_vals = {10.0};
  std::vector<f64> const alt_vals = {1.0, 2.0, 3.0};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE(result.has_value());
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  CHECK(result.value() < 0.0);
}

TEST_CASE("MannWhitneyEffectSize handles N=2 minimum (1 ref, 1 alt)",
          "[lancet][base][MannWhitneyEffectSize]") {
  // Smallest legal sample. With N=2 the variance denominator N*(N-1)=2 stays
  // positive only because tie_correction is 0; if ties were present here the
  // variance would underflow. Verify a value-distinct pair returns a finite
  // engaged optional.
  std::vector<f64> const ref_vals = {1.0};
  std::vector<f64> const alt_vals = {2.0};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE(result.has_value());
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  CHECK(std::isfinite(result.value()));
}

TEST_CASE("MannWhitneyEffectSize handles asymmetric group sizes",
          "[lancet][base][MannWhitneyEffectSize]") {
  // 1 alt vs 9 ref — the kind of imbalance somatic callers see at low VAF.
  std::vector<f64> const ref_vals = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  std::vector<f64> const alt_vals = {20.0};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE(result.has_value());
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  CHECK(result.value() > 0.0);
}

TEST_CASE("MannWhitneyEffectSize handles 1000-element groups (max-size sanity)",
          "[lancet][base][MannWhitneyEffectSize]") {
  // Exercises the rank-sum and tie-correction loops at scale. ALT sample is
  // shifted strictly above REF, so the effect size must be positive.
  std::vector<f64> ref_vals;
  std::vector<f64> alt_vals;
  ref_vals.reserve(1000);
  alt_vals.reserve(1000);
  for (i32 idx = 0; idx < 1000; ++idx) {
    ref_vals.push_back(static_cast<f64>(idx));
    alt_vals.push_back(static_cast<f64>(idx + 1000));
  }
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE(result.has_value());
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  CHECK(result.value() > 0.0);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  Permutation invariance                                                  ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("MannWhitneyEffectSize is invariant under permutation of ref or alt within their groups",
          "[lancet][base][MannWhitneyEffectSize]") {
  // Mann-Whitney is rank-based, so the order of values inside a group has no
  // effect on the rank sum. A regression that accidentally indexed by position
  // rather than value would surface here.
  std::vector<f64> const ref_sorted = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<f64> const ref_shuffled = {3.0, 5.0, 1.0, 4.0, 2.0};
  std::vector<f64> const alt_sorted = {6.0, 7.0, 8.0, 9.0, 10.0};
  std::vector<f64> const alt_shuffled = {9.0, 6.0, 8.0, 10.0, 7.0};

  auto const baseline =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_sorted), absl::MakeConstSpan(alt_sorted));
  auto const ref_permuted = MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_shuffled),
                                                       absl::MakeConstSpan(alt_sorted));
  auto const alt_permuted = MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_sorted),
                                                       absl::MakeConstSpan(alt_shuffled));
  auto const both_permuted = MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_shuffled),
                                                        absl::MakeConstSpan(alt_shuffled));

  REQUIRE(baseline.has_value());
  REQUIRE(ref_permuted.has_value());
  REQUIRE(alt_permuted.has_value());
  REQUIRE(both_permuted.has_value());
  // engaged-optional asserted on the prior REQUIREs; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  CHECK(ref_permuted.value() == Catch::Approx(baseline.value()).margin(1e-12));
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  CHECK(alt_permuted.value() == Catch::Approx(baseline.value()).margin(1e-12));
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  CHECK(both_permuted.value() == Catch::Approx(baseline.value()).margin(1e-12));
}

TEST_CASE("MannWhitneyEffectSize swapping ref<->alt negates the effect size",
          "[lancet][base][MannWhitneyEffectSize]") {
  // The U statistic is asymmetric: U_ref + U_alt = m·n. Dividing by √Var(U)
  // and √N preserves this antisymmetry — swapping the inputs negates the
  // returned effect size. A regression that inadvertently took |Z| would
  // surface here.
  std::vector<f64> const low_vals = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<f64> const high_vals = {6.0, 7.0, 8.0, 9.0, 10.0};
  auto const positive =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(low_vals), absl::MakeConstSpan(high_vals));
  auto const negative =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(high_vals), absl::MakeConstSpan(low_vals));

  REQUIRE(positive.has_value());
  REQUIRE(negative.has_value());
  // engaged-optional asserted on the prior REQUIREs; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  CHECK(positive.value() == Catch::Approx(-negative.value()).margin(1e-12));
}

}  // namespace lancet::base::tests
