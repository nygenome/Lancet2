#ifndef SRC_LANCET_BASE_MANN_WHITNEY_H_
#define SRC_LANCET_BASE_MANN_WHITNEY_H_

// ============================================================================
// Mann-Whitney U Test — Effect Size (Cohen's d analog)
//
// Provides a non-parametric effect size measure for whether two independent
// samples are drawn from the same distribution, derived from the standard
// Mann-Whitney U test (Wilcoxon Rank-Sum Test).
//
// Unlike the raw Z-score (which scales with √N due to the Central Limit
// Theorem), the effect size Z/√N is coverage-invariant: the same biological
// bias produces the same value regardless of sequencing depth. This makes
// it suitable for ML models that must generalize across coverage regimes.
//
// ── Mathematical Foundation ─────────────────────────────────────────────────
//
// Given two samples X = {x_1, ..., x_m} (REF) and Y = {y_1, ..., y_n} (ALT):
//
//   1. Pool all N = m + n observations and assign ranks 1..N (ascending).
//      Ties receive the mean rank of the tied group (mid-rank method).
//
//   2. Compute the U statistic for the ALT sample:
//
//        U_alt = R_alt - n*(n+1)/2
//
//      where R_alt is the sum of ranks assigned to ALT observations.
//
//   3. Under H₀ (both samples from the same distribution):
//
//        E[U] = m * n / 2
//
//   4. Variance with tie correction (Lehmann, 2006):
//
//        Var(U) = (m * n / 12) * [(N + 1) - Σₖ(tₖ³ - tₖ) / (N * (N - 1))]
//
//   5. Raw Z-score: Z = (U - E[U]) / √Var(U)
//
//   6. Effect size normalization: Z / √N where N = m + n.
//
//      Since Z ∝ √N × δ for a constant effect size δ (by CLT),
//      dividing by √N recovers δ — the standardized effect size
//      (analogous to Cohen's d for parametric tests).
//
// ── Coverage Stability ──────────────────────────────────────────────────────
//
// For a constant mild bias (e.g., ALT MAPQ 2 units lower than REF):
//
//   20×:   raw Z ≈ −1.5,  Z/√N ≈ −0.34
//   60×:   raw Z ≈ −2.6,  Z/√N ≈ −0.34
//   100×:  raw Z ≈ −3.3,  Z/√N ≈ −0.33
//   1000×: raw Z ≈ −10.5, Z/√N ≈ −0.33
//   2000×: raw Z ≈ −14.9, Z/√N ≈ −0.33
//
// The effect size varies < 3% across 20×–2000× for any fixed bias level.
//
// ── Edge Cases ──────────────────────────────────────────────────────────────
//
// When the test cannot be computed (one or both groups empty, or all values
// identical producing zero variance), 0.0 is returned. This is correct:
// with no observable difference between groups, there is no measured bias.
//
// Common biological scenarios producing 0.0:
//   - Homozygous ALT (1/1): all reads support ALT, 0 REF observations
//   - Allelic dropout: all reads randomly sample one allele
//   - Homozygous REF (0/0): 0 ALT observations (no variant reads)
//   - Identical values: all reads have the same MAPQ/BQ/position
//
// ── References ──────────────────────────────────────────────────────────────
//
//   - Mann, H.B. & Whitney, D.R. (1947). Annals of Math. Statistics, 18(1).
//   - Lehmann, E.L. (2006). Nonparametrics: Statistical Methods Based on Ranks.
//   - Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences.
//
// ============================================================================

#include "lancet/base/types.h"

#include "absl/types/span.h"

#include <algorithm>
#include <vector>

#include <cmath>

namespace lancet::base {

// Computes the Mann-Whitney U effect size (Z/√N) comparing two independent
// samples. This is the coverage-normalized analog of the Z-score.
//
// Returns positive values if alt_vals tend to be higher, negative if lower.
// Returns 0.0 when either sample is empty (test not applicable) or when
// all values are identical (zero variance, no discriminating power).
//
// Template parameter T must be an arithmetic type (u8, i32, f64, etc.).
// Values are promoted to f64 internally for rank computation.
template <typename T>
[[nodiscard]] auto MannWhitneyEffectSize(absl::Span<T const> ref_vals, absl::Span<T const> alt_vals)
    -> f64 {
  // Require at least one observation in each group. Without both groups,
  // there is no observable bias — return 0.0.
  if (ref_vals.empty() || alt_vals.empty()) {
    return 0.0;
  }

  auto const n_ref = static_cast<f64>(ref_vals.size());
  auto const n_alt = static_cast<f64>(alt_vals.size());
  auto const total = ref_vals.size() + alt_vals.size();

  // ── Step 1: Pool and sort ────────────────────────────────────────────────
  // Tag each value with its group membership for rank assignment.
  struct TaggedValue {
    f64 mValue;
    bool mIsAlt;
  };

  std::vector<TaggedValue> pooled;
  pooled.reserve(total);
  for (auto const val : ref_vals) {
    pooled.push_back({static_cast<f64>(val), false});
  }
  for (auto const val : alt_vals) {
    pooled.push_back({static_cast<f64>(val), true});
  }

  std::sort(pooled.begin(), pooled.end(),
            [](TaggedValue const& lhs, TaggedValue const& rhs) -> bool {
              return lhs.mValue < rhs.mValue;
            });

  // ── Step 2: Assign mid-ranks and accumulate ALT rank sum + tie term ──────
  // Tied values receive the average of the ranks they would span.
  // Example: values [3, 5, 5, 7] get ranks [1, 2.5, 2.5, 4].
  f64 alt_rank_sum = 0.0;
  f64 tie_correction = 0.0;  // Σ(tₖ³ - tₖ) accumulated over all tie groups

  for (usize i = 0; i < total;) {
    // Find the extent of this tie group: [i, jdx)
    usize jdx = i;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
    while (jdx < total && pooled[jdx].mValue == pooled[i].mValue) {
      ++jdx;
    }

    // Mid-rank for positions [i+1, jdx] (1-indexed): (i+1 + jdx) / 2
    auto const mid_rank = static_cast<f64>(i + 1 + jdx) / 2.0;

    // Tie group size for variance correction
    auto const tie_size = static_cast<f64>(jdx - i);
    tie_correction += ((tie_size * tie_size * tie_size) - tie_size);

    // Accumulate rank sum for ALT observations only
    for (usize k = i; k < jdx; ++k) {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
      if (pooled[k].mIsAlt) {
        alt_rank_sum += mid_rank;
      }
    }

    i = jdx;
  }

  // ── Step 3: Compute U statistic ──────────────────────────────────────────
  // U_alt = R_alt - n_alt*(n_alt+1)/2
  // Counts the number of (ref, alt) pairs where ref value > alt value.
  auto const u_stat = alt_rank_sum - ((n_alt * (n_alt + 1.0)) / 2.0);

  // ── Step 4: Expected value and variance under H₀ ─────────────────────────
  auto const mean_u = (n_ref * n_alt) / 2.0;
  auto const n_total = static_cast<f64>(total);

  // Variance with tie correction (Lehmann, 2006, §1.4):
  //   Var(U) = (m*n/12) * [(N+1) - Σ(t³-t) / (N*(N-1))]
  // The tie term reduces variance because tied values carry less information.
  auto const var_u =
      (n_ref * n_alt / 12.0) * ((n_total + 1.0) - (tie_correction / (n_total * (n_total - 1.0))));

  // Zero variance means all values are identical — no observable bias.
  if (var_u <= 0.0) {
    return 0.0;
  }

  // ── Step 5: Effect size = Z / √N ─────────────────────────────────────────
  // Raw Z-score divided by √N to remove √N power amplification.
  // This recovers the standardized effect size (Cohen's d analog).
  auto const z_score = (u_stat - mean_u) / std::sqrt(var_u);
  return z_score / std::sqrt(n_total);
}

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_MANN_WHITNEY_H_
