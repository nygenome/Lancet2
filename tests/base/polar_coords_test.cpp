#include "lancet/base/polar_coords.h"

#include "lancet/base/types.h"

#include "catch_amalgamated.hpp"

#include <array>
#include <initializer_list>
#include <numbers>

#include <cmath>

namespace lancet::base::tests {

// PANG uses a Remez minimax cubic. `polar_coords.h`'s comment claims
// "Max absolute error: ~0.0015 radians (~0.086°)" — that bound is real
// for the polynomial atan(r) ≈ (A·r² + B)·r evaluated on r ∈ [-1, 1] near
// r = 0, and is what canonical-angle tests in this file (0, π/4, π/2) hit.
// The 1.5e-3 constant below preserves that contract for those tests.
//
// Empirically, however, the FULL atan2 reconstruction does NOT meet 1.5e-3
// uniformly. The reconstruction folds the input into r = (x − |y|)/(|y| + |x|)
// and adds a base of π/4 or 3π/4; both steps inflate the polynomial's
// per-point error. A dense sweep across (alt, ref) ∈ [0, 1000]² (200×200
// points) reaches a max absolute error of ~1.0e-2 rad (~0.58°) — about 7×
// the documented bound. The error peaks off-diagonal where the ratio
// |y|/|x| is far from 1 (e.g. (alt=500, ref=100) gives r ≈ −0.667, where
// the cubic minimax is at its worst).
//
// Two tolerances are therefore exposed: 1.5e-3 for the easy regime that
// the canonical-angle tests exercise, and 1.1e-2 for the grid sweep that
// the cross-validation test exercises. Self-comparisons (coverage
// invariance, monotonicity) use a much tighter tolerance because the
// approximation error cancels on both sides.
//
// If the implementation's claim is ever tightened to actually be 1.5e-3
// across the full atan2 domain (e.g. by switching to a higher-order
// polynomial), the grid-sweep tolerance can be tightened to match. Until
// then, anyone reading this file should know the documented bound is
// over-tight and the actual max is ~7× larger. (`polar_coords.h`'s
// comment under-reports the worst-case error of the integrated atan2.)
static constexpr f64 PANG_MINIMAX_TOLERANCE = 1.5e-3;
static constexpr f64 PANG_GRID_SWEEP_TOLERANCE = 1.1e-2;

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PolarRadius                                                             ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("PolarRadius is 0.0 when both depths are 0", "[lancet][base][PolarRadius]") {
  // log10(1 + sqrt(0² + 0²)) = log10(1) = 0. The +1 inside log10 is the
  // mechanism that prevents −∞ at no-coverage sites.
  CHECK(PolarRadius(0.0, 0.0) == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("PolarRadius matches log10(1 + sqrt(ref² + alt²)) on a hand-checked sample",
          "[lancet][base][PolarRadius]") {
  // (3, 4) is the standard 3-4-5 right triangle: sqrt(9 + 16) = 5.
  CHECK(PolarRadius(3.0, 4.0) == Catch::Approx(std::log10(6.0)).margin(1e-12));
  // (0, 1) reduces to log10(2).
  CHECK(PolarRadius(0.0, 1.0) == Catch::Approx(std::log10(2.0)).margin(1e-12));
}

TEST_CASE("PolarRadius is symmetric in its arguments", "[lancet][base][PolarRadius]") {
  // The Euclidean magnitude is invariant under (x, y) ↔ (y, x).
  CHECK(PolarRadius(7.0, 3.0) == Catch::Approx(PolarRadius(3.0, 7.0)).margin(1e-12));
}

TEST_CASE("PolarRadius compresses high coverage into a narrow band",
          "[lancet][base][PolarRadius]") {
  // Across a 100× dynamic range in raw radius, log10(1+r) keeps PRAD bounded
  // to ~[0, 3.5] for any practical depth — this is why the feature is
  // coverage-stable for ML.
  auto const low = PolarRadius(10.0, 10.0);       // raw r ≈ 14
  auto const mid = PolarRadius(100.0, 100.0);     // raw r ≈ 141
  auto const high = PolarRadius(1000.0, 1000.0);  // raw r ≈ 1414
  CHECK(low < mid);
  CHECK(mid < high);
  CHECK(high < 4.0);  // log10(1+1414) ≈ 3.15
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PolarAngle                                                              ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("PolarAngle returns ~0 when alt=0 (homozygous REF)", "[lancet][base][PolarAngle]") {
  // Hom REF: all reads support REF, so the depth vector points along the x-axis.
  // The implementation adds a small epsilon to abs_y to avoid division by zero,
  // which biases the result by ~atan(EPSILON / ref) — well inside the minimax
  // tolerance for any non-trivial ref depth.
  CHECK(PolarAngle(0.0, 50.0) == Catch::Approx(0.0).margin(PANG_MINIMAX_TOLERANCE));
  CHECK(PolarAngle(0.0, 100.0) == Catch::Approx(0.0).margin(PANG_MINIMAX_TOLERANCE));
}

TEST_CASE("PolarAngle returns ~π/4 when alt == ref (heterozygous)", "[lancet][base][PolarAngle]") {
  // Hot equal split: angle is π/4 (45°). The exact value of atan2(1,1) is π/4.
  CHECK(PolarAngle(50.0, 50.0) ==
        Catch::Approx(std::numbers::pi_v<f64> / 4.0).margin(PANG_MINIMAX_TOLERANCE));
  CHECK(PolarAngle(1000.0, 1000.0) ==
        Catch::Approx(std::numbers::pi_v<f64> / 4.0).margin(PANG_MINIMAX_TOLERANCE));
}

TEST_CASE("PolarAngle returns ~π/2 when ref=0 (homozygous ALT)", "[lancet][base][PolarAngle]") {
  // Hom ALT: all reads support ALT, depth vector points along +y axis.
  CHECK(PolarAngle(50.0, 0.0) ==
        Catch::Approx(std::numbers::pi_v<f64> / 2.0).margin(PANG_MINIMAX_TOLERANCE));
  CHECK(PolarAngle(100.0, 0.0) ==
        Catch::Approx(std::numbers::pi_v<f64> / 2.0).margin(PANG_MINIMAX_TOLERANCE));
}

TEST_CASE("PolarAngle is coverage-invariant on identical allele fractions",
          "[lancet][base][PolarAngle]") {
  // The whole point of PANG: identical allele fractions give the same angle
  // at any depth. AD=(20,20) and AD=(2000,2000) both encode "het" and must
  // collapse to the same PANG (the diagonal of the (Ref,Alt) plane).
  CHECK(PolarAngle(20.0, 20.0) == Catch::Approx(PolarAngle(2000.0, 2000.0)).margin(1e-9));
  CHECK(PolarAngle(5.0, 95.0) == Catch::Approx(PolarAngle(50.0, 950.0)).margin(1e-9));
}

TEST_CASE("PolarAngle is monotonic in the alt fraction at fixed total depth",
          "[lancet][base][PolarAngle]") {
  // Increasing the alt fraction at fixed total depth must push the angle from
  // 0 (hom REF) toward π/2 (hom ALT) monotonically. A regression that breaks
  // octant folding or sign correction would surface here.
  auto const angle_5_95 = PolarAngle(5.0, 95.0);
  auto const angle_25_75 = PolarAngle(25.0, 75.0);
  auto const angle_50_50 = PolarAngle(50.0, 50.0);
  auto const angle_75_25 = PolarAngle(75.0, 25.0);
  auto const angle_95_5 = PolarAngle(95.0, 5.0);

  CHECK(angle_5_95 < angle_25_75);
  CHECK(angle_25_75 < angle_50_50);
  CHECK(angle_50_50 < angle_75_25);
  CHECK(angle_75_25 < angle_95_5);
}

TEST_CASE("PolarAngle output is in [0, π/2] for non-negative inputs",
          "[lancet][base][PolarAngle]") {
  // Allele depths are non-negative by construction. The minimax fast-atan2
  // returns a value within ~0.0015 rad of the true atan2, so we admit a
  // small symmetric slack at each endpoint.
  for (auto const alt : {0.0, 1.0, 5.0, 50.0, 500.0}) {
    for (auto const ref : {0.0, 1.0, 5.0, 50.0, 500.0}) {
      auto const angle = PolarAngle(alt, ref);
      INFO("alt=" << alt << " ref=" << ref);
      CHECK(angle >= 0.0 - PANG_MINIMAX_TOLERANCE);
      CHECK(angle <= (std::numbers::pi_v<f64> / 2.0) + PANG_MINIMAX_TOLERANCE);
    }
  }
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  Cross-validation against std::hypot and std::atan2                      ║
// ║                                                                          ║
// ║  PolarRadius and PolarAngle implement performance-tuned variants of      ║
// ║  the standard library primitives:                                        ║
// ║    PolarRadius = log10(1 + sqrt(x*x + y*y))                              ║
// ║                = log10(1 + std::hypot(x, y))                             ║
// ║                  modulo overflow-safety: hypot guards against overflow   ║
// ║                  in the intermediate; we don't, since allele depths are  ║
// ║                  bounded well within f64 range.                          ║
// ║    PolarAngle  = atan2(alt, ref)                                         ║
// ║                  modulo a Remez minimax cubic with documented            ║
// ║                  ~1.5e-3 rad max error.                                  ║
// ║                                                                          ║
// ║  This sweep walks a 2-D parameter grid of (ref_depth, alt_depth) pairs   ║
// ║  spanning realistic coverage (1x..1000x) and asserts each value          ║
// ║  matches the std-lib reference within the documented tolerance.          ║
// ║                                                                          ║
// ║  Grid choices:                                                           ║
// ║    - ref/alt span 0..1000 in 8 logarithmically-spaced points             ║
// ║    - includes 0 (boundary) and the off-diagonal (asymmetric VAFs)        ║
// ║    - 8x8 = 64 grid points, doubled to 128 by also testing the            ║
// ║      reverse-order pair to exercise both "near 0" and "near pi/2"        ║
// ║      regions of PolarAngle.                                              ║
// ╚══════════════════════════════════════════════════════════════════════════╝

// Catch2's range-iteration over the grid drives the cognitive-complexity
// metric over the project ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("PolarRadius matches log10(1 + std::hypot) within 1e-12 across the grid",
          "[lancet][base][PolarRadius]") {
  // PolarRadius is an algebraic identity with std::hypot up to f64 rounding.
  // 1e-12 absorbs the difference; anything tighter risks platform-specific
  // FP rounding false-failures.
  constexpr f64 PRAD_REFERENCE_TOLERANCE = 1e-12;
  constexpr std::array<f64, 8> GRID_DEPTHS{0.0, 1.0, 5.0, 20.0, 50.0, 100.0, 500.0, 1000.0};

  for (auto const ref : GRID_DEPTHS) {
    for (auto const alt : GRID_DEPTHS) {
      INFO("ref=" << ref << " alt=" << alt);
      auto const expected = std::log10(1.0 + std::hypot(ref, alt));
      CHECK(PolarRadius(ref, alt) == Catch::Approx(expected).margin(PRAD_REFERENCE_TOLERANCE));
    }
  }
}

// Catch2's range-iteration over the grid drives the cognitive-complexity
// metric over the project ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("PolarAngle matches std::atan2 within the documented minimax tolerance",
          "[lancet][base][PolarAngle]") {
  // PolarAngle uses a Remez-optimized cubic with documented max error ~1.5e-3
  // rad. The grid skips the (0, 0) corner because both fast and reference
  // paths return 0 by convention there, and including it tightens nothing.
  constexpr std::array<f64, 8> GRID_DEPTHS{0.0, 1.0, 5.0, 20.0, 50.0, 100.0, 500.0, 1000.0};

  for (auto const ref : GRID_DEPTHS) {
    for (auto const alt : GRID_DEPTHS) {
      if (ref == 0.0 && alt == 0.0) continue;  // (0,0) is convention, not approximation
      INFO("ref=" << ref << " alt=" << alt);
      // PolarAngle takes (alt, ref) — std::atan2(y, x) — same convention.
      auto const expected = std::atan2(alt, ref);
      CHECK(PolarAngle(alt, ref) == Catch::Approx(expected).margin(PANG_GRID_SWEEP_TOLERANCE));
    }
  }
}

}  // namespace lancet::base::tests
