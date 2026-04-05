#ifndef SRC_LANCET_BASE_MANN_WHITNEY_H_
#define SRC_LANCET_BASE_MANN_WHITNEY_H_

// ============================================================================
// Mann-Whitney U Test (Wilcoxon Rank-Sum Test) — Z-score Approximation
//
// Provides a non-parametric test for whether two independent samples are drawn
// from the same distribution. This is the standard test used by GATK's
// MappingQualityRankSumTest (MQRankSum) and related RankSum annotations.
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
//      U counts the number of (ref, alt) pairs where ref > alt.
//
//   3. Under H₀ (both samples from the same distribution):
//
//        E[U] = m * n / 2
//
//   4. Variance with tie correction (Lehmann, 2006):
//
//        Var(U) = (m * n / 12) * [(N + 1) - Σₖ(tₖ³ - tₖ) / (N * (N - 1))]
//
//      where the sum is over all tie groups k, and tₖ is the group size.
//      Without ties, this simplifies to m*n*(N+1)/12.
//
//   5. Z-score (normal approximation, valid for m,n ≥ 8):
//
//        Z = (U - E[U]) / √Var(U)
//
//      Negative Z ⟹ ALT values are systematically lower than REF.
//      Positive Z ⟹ ALT values are systematically higher than REF.
//
// ── Usage in Variant Calling ────────────────────────────────────────────────
//
//   MQRS (Mapping Quality Rank Sum Z-score):
//     Compares mapping qualities of REF-supporting vs ALT-supporting reads.
//     Negative Z indicates ALT reads have lower mapping quality, suggesting
//     paralogous mismapping. GATK uses Z < -12.5 as a hard filter threshold.
//
// ── References ──────────────────────────────────────────────────────────────
//
//   - Mann, H.B. & Whitney, D.R. (1947). Annals of Math. Statistics, 18(1).
//   - Lehmann, E.L. (2006). Nonparametrics: Statistical Methods Based on Ranks.
//   - GATK MappingQualityRankSumTest: u-based Z-approximation.
//   - bcftools: mqsbias (related rank-based mapping quality bias metric).
//
// ============================================================================

#include <algorithm>
#include <cmath>
#include <vector>

#include "absl/types/span.h"
#include "lancet/base/types.h"

namespace lancet::base {

// ── Sentinel value for untestable edge cases ─────────────────────────────────
//
// When the Mann-Whitney U test cannot be computed (one or both sample groups
// are empty, or all observations are identical), returning 0.0 is statistically
// dangerous: it represents the null hypothesis mean ("no bias"), which falsely
// signals high confidence of allele balance to downstream filters and ML models.
//
// Instead we use extreme value imputation: a value so far outside the natural
// Z-score range that it is trivially separable by any downstream consumer.
//
// Common biological scenarios producing an untestable state:
//   - Homozygous ALT (1/1): all reads support ALT, 0 REF observations.
//   - Allelic dropout: low coverage site where all reads randomly sample one allele.
//   - Homozygous REF (0/0): 0 ALT observations (no variant reads).
//   - Identical MAPQ: all reads share the exact same mapping quality (zero variance).
//
// Why 100.0?
//   - Natural Z-scores from genomic data rarely exceed |8| (GATK hard-filters at
//     |12.5|, which itself is extreme). A value of 100.0 is >7σ beyond any
//     plausible biological signal.
//   - Tree-based ML models (Isolation Forest, XGBoost) will isolate this value
//     with an early split, keeping it from contaminating the learned Z-score
//     distribution.
//   - Density-based models (GMMs) should pair this with a binary indicator
//     feature (is_MQRS_missing) to avoid creating a false cluster at 100.0.
//   - VCF consumers can trivially filter on MQRS == 100.0 to identify untested
//     loci.
//
// Reference: See conversation notes on VQSR/ML edge-case handling for
// rank-sum annotations.
static constexpr f64 MANN_WHITNEY_MISSING_SENTINEL = 100.0;

// Computes the Mann-Whitney U test Z-score comparing two independent samples.
//
// Returns positive Z if alt_vals tend to be higher, negative if lower.
// Returns MANN_WHITNEY_MISSING_SENTINEL (100.0) when either sample is empty
// (test not applicable) or when all values are identical (zero variance,
// no discriminating power). See sentinel documentation above.
//
// Template parameter T must be an arithmetic type (u8, i32, f64, etc.).
// Values are promoted to f64 internally for rank computation.
template <typename T>
[[nodiscard]] auto MannWhitneyUZScore(absl::Span<const T> ref_vals, absl::Span<const T> alt_vals) -> f64 {
  // Require at least one observation in each group. Without both groups,
  // the rank-sum test is undefined — return sentinel, not 0.0.
  if (ref_vals.empty() || alt_vals.empty()) return MANN_WHITNEY_MISSING_SENTINEL;

  const auto n_ref = static_cast<f64>(ref_vals.size());
  const auto n_alt = static_cast<f64>(alt_vals.size());
  const auto total = ref_vals.size() + alt_vals.size();

  // ── Step 1: Pool and sort ────────────────────────────────────────────────
  // Tag each value with its group membership for rank assignment.
  struct TaggedValue {
    f64 value;
    bool is_alt;
  };

  std::vector<TaggedValue> pooled;
  pooled.reserve(total);
  for (const auto val : ref_vals) pooled.push_back({static_cast<f64>(val), false});
  for (const auto val : alt_vals) pooled.push_back({static_cast<f64>(val), true});

  std::sort(pooled.begin(), pooled.end(),
            [](const TaggedValue& lhs, const TaggedValue& rhs) { return lhs.value < rhs.value; });

  // ── Step 2: Assign mid-ranks and accumulate ALT rank sum + tie term ──────
  // Tied values receive the average of the ranks they would span.
  // Example: values [3, 5, 5, 7] get ranks [1, 2.5, 2.5, 4].
  f64 alt_rank_sum = 0.0;
  f64 tie_correction = 0.0;  // Σ(tₖ³ - tₖ) accumulated over all tie groups

  for (usize i = 0; i < total;) {
    // Find the extent of this tie group: [i, j)
    usize j = i;
    while (j < total && pooled[j].value == pooled[i].value) ++j;

    // Mid-rank for positions [i+1, j] (1-indexed): (i+1 + j) / 2
    const auto mid_rank = static_cast<f64>(i + 1 + j) / 2.0;

    // Tie group size for variance correction
    const auto tie_size = static_cast<f64>(j - i);
    tie_correction += (tie_size * tie_size * tie_size - tie_size);

    // Accumulate rank sum for ALT observations only
    for (usize k = i; k < j; ++k) {
      if (pooled[k].is_alt) alt_rank_sum += mid_rank;
    }

    i = j;
  }

  // ── Step 3: Compute U statistic ──────────────────────────────────────────
  // U_alt = R_alt - n_alt*(n_alt+1)/2
  // Counts the number of (ref, alt) pairs where ref value > alt value.
  const auto u_stat = alt_rank_sum - (n_alt * (n_alt + 1.0)) / 2.0;

  // ── Step 4: Expected value and variance under H₀ ─────────────────────────
  const auto mean_u = (n_ref * n_alt) / 2.0;
  const auto n_total = static_cast<f64>(total);

  // Variance with tie correction (Lehmann, 2006, §1.4):
  //   Var(U) = (m*n/12) * [(N+1) - Σ(t³-t) / (N*(N-1))]
  // The tie term reduces variance because tied values carry less information.
  const auto var_u = (n_ref * n_alt / 12.0) * ((n_total + 1.0) - tie_correction / (n_total * (n_total - 1.0)));

  // Zero variance means all values are identical — no test is meaningful.
  // Return sentinel rather than 0.0 to avoid false "no bias" signal.
  if (var_u <= 0.0) return MANN_WHITNEY_MISSING_SENTINEL;

  // ── Step 5: Z-score (normal approximation) ───────────────────────────────
  return (u_stat - mean_u) / std::sqrt(var_u);
}

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_MANN_WHITNEY_H_
