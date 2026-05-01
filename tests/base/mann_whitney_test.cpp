#include "lancet/base/mann_whitney.h"

#include "absl/types/span.h"
#include "catch_amalgamated.hpp"

#include <vector>

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
          "[lancet][base][MannWhitney]") {
  std::vector<f64> const ref_vals = {};
  std::vector<f64> const alt_vals = {1.0, 2.0, 3.0};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE_FALSE(result.has_value());
}

TEST_CASE("MannWhitneyEffectSize returns nullopt when ALT group is empty",
          "[lancet][base][MannWhitney]") {
  std::vector<f64> const ref_vals = {1.0, 2.0, 3.0};
  std::vector<f64> const alt_vals = {};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE_FALSE(result.has_value());
}

TEST_CASE("MannWhitneyEffectSize returns nullopt when both groups are empty",
          "[lancet][base][MannWhitney]") {
  std::vector<f64> const ref_vals = {};
  std::vector<f64> const alt_vals = {};
  auto const result =
      MannWhitneyEffectSize<f64>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE_FALSE(result.has_value());
}

TEST_CASE("MannWhitneyEffectSize returns 0.0 when all values are identical",
          "[lancet][base][MannWhitney]") {
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
          "[lancet][base][MannWhitney]") {
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

TEST_CASE("MannWhitneyEffectSize works with u8 type", "[lancet][base][MannWhitney]") {
  std::vector<u8> const ref_vals = {60, 60, 60};
  std::vector<u8> const alt_vals = {};
  auto const result =
      MannWhitneyEffectSize<u8>(absl::MakeConstSpan(ref_vals), absl::MakeConstSpan(alt_vals));
  REQUIRE_FALSE(result.has_value());
}

}  // namespace lancet::base::tests
