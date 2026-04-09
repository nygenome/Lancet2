#ifndef SRC_LANCET_CBDG_GRAPH_COMPLEXITY_H_
#define SRC_LANCET_CBDG_GRAPH_COMPLEXITY_H_

#include "lancet/base/types.h"

#include <cmath>

namespace lancet::cbdg {

class Graph;  // forward declaration — avoids circular dependency with graph.h

// ============================================================================
// GraphComplexity — O(V+E) topology metrics for a single graph component
//
// Characterizes the structural complexity of a local assembly graph component.
// Used for two purposes:
//
//   1. Walk enumeration gating: IsComplex() triggers k-mer upsizing when the
//      graph is too tangled for efficient path enumeration.
//
//   2. ML feature generation: Metrics are compressed into the Graph
//      Entanglement Index (GEI) and two orthogonal features for downstream
//      variant classification. See GraphEntanglementIndex() method below.
//
// ── Retained Metrics ────────────────────────────────────────────────────────
//
//  CyclomaticComplexity (CC = E - V + 1):
//    Number of independent cycles. CC=0 → linear chain, CC=1 → single
//    variant bubble, CC≥5 → tangled repeat region. Used in both IsComplex()
//    and as a GEI numerator term.
//
//  NumBranchPoints:
//    Nodes with ≥2 edges in some direction. 2 per clean variant (entry
//    and exit of a bubble), ≥5 → anomalous branching. Used in IsComplex()
//    and as a GEI numerator term.
//
//  UnitigRatio:
//    Fraction of nodes with exactly 1-in/1-out (linear chain nodes).
//    >90% = high quality, <70% = complex topology. GEI denominator term.
//
//  CoverageCv (σ/μ of TotalReadSupport):
//    Coefficient of variation of total read support across nodes.
//    Low CV = uniform coverage, high CV = collapsed or divergent region.
//    GEI numerator term — distinguishes true polymorphism from collapsed
//    repeats.
//
//  MaxSingleDirDegree:
//    Maximum outgoing edges in any single sign direction. ≤3 normal,
//    >5 = hub k-mer causing BFS blowup. Orthogonal VCF feature.
//
//  TipToPathCovRatio:
//    Mean dead-end node coverage / mean linear (unitig) node coverage.
//    ≈0 = clean (tips have negligible coverage), high = assembly tearing.
//    Orthogonal VCF feature.
//
// Quality interpretation:
//   HIGH:   CC ≤ 2, BP ≤ 4, UnitigRatio > 0.90, CovCV < 0.5
//   MEDIUM: CC 3–10, BP 5–15, UnitigRatio 0.70–0.90, CovCV 0.5–1.5
//   LOW:    CC > 10, BP > 15, UnitigRatio < 0.70, CovCV > 1.5
//
//  NumNodes / NumEdges:
//    Used as intermediates during CC computation — not stored on the class.
//
// ============================================================================
class GraphComplexity {
 public:
  // ── Accessors ─────────────────────────────────────────────────────────
  [[nodiscard]] auto CyclomaticComplexity() const noexcept -> usize {
    return mCyclomaticComplexity;
  }
  [[nodiscard]] auto NumBranchPoints() const noexcept -> usize { return mNumBranchPoints; }
  [[nodiscard]] auto UnitigRatio() const noexcept -> f64 { return mUnitigRatio; }
  [[nodiscard]] auto CoverageCv() const noexcept -> f64 { return mCoverageCv; }
  [[nodiscard]] auto MaxSingleDirDegree() const noexcept -> usize { return mMaxSingleDirDegree; }
  [[nodiscard]] auto TipToPathCovRatio() const noexcept -> f64 { return mTipToPathCovRatio; }

  /// Thresholds for classifying a graph as too complex for efficient walk
  /// enumeration at the current k-mer size. When IsComplex() returns true,
  /// the caller retries with a larger k to collapse branches.
  ///
  /// Biological interpretation (in a ~1000bp window with k=25, ~975 nodes):
  ///   CC = 50 ≈ 50 independent variant signals in 1000bp (~1 per 20bp),
  ///   which is ~50× the normal human heterozygosity rate (~1 SNP/1000bp).
  ///   Br = 50 ≈ 5% of all positions forking into multiple continuations.
  ///   Together, CC>=50 AND Br>=50 indicates the graph is not driven by
  ///   clean biological variation but by tandem repeats (STRs/VNTRs),
  ///   segmental duplications, structural variant breakpoints, or mapping
  ///   artifacts. A larger k collapses short repeats and merges paralogous
  ///   paths, simplifying the graph toward its true biological structure.
  ///
  /// Derived from whole-chromosome profiling (chr4, 233,479 windows, 235,533
  /// component entries).
  ///
  ///   Threshold analysis (AND — both conditions must hold):
  ///     CC>=50  AND Br>=50:   146 wins, avg  5.8s | 233,333 normal, avg 414ms (14×)
  ///     CC>=50  AND Br>=100:   76 wins, avg 10.4s | 233,403 normal, avg 414ms (25×)
  ///     CC>=100 AND Br>=100:   31 wins, avg 15.2s | 233,448 normal, avg 415ms (37×)
  ///     CC>=100 AND Br>=200:   18 wins, avg 25.0s | 233,461 normal, avg 415ms (60×)
  ///
  ///   Using AND (vs OR) avoids skipping windows where only one metric is
  ///   elevated — e.g. a high-CC but low-branch graph may still enumerate
  ///   quickly. Both metrics must agree that the graph is pathological.
  ///
  ///   Distribution stats (235,533 components):
  ///     CC:  median=3, p95=9, p99=16, max=487
  ///     Br:  median=6, p95=32, p99=32, max=1127
  ///
  static constexpr usize MAX_CYCLOMATIC_COMPLEXITY = 50;
  static constexpr usize MAX_BRANCH_POINTS = 50;

  /// Returns true if the graph is too complex for efficient walk enumeration.
  /// When true, the caller should retry with a larger k-mer size to collapse
  /// branches and reduce cyclomatic complexity.
  [[nodiscard]] constexpr auto IsComplex() const -> bool {
    return mCyclomaticComplexity >= MAX_CYCLOMATIC_COMPLEXITY &&
           mNumBranchPoints >= MAX_BRANCH_POINTS;
  }

  // ── Graph Entanglement Index (GEI) ──────────────────────────────────────
  //
  // Compresses four collinear topology metrics (CC, BP, CoverageCv,
  // UnitigRatio) into a single continuous scalar for ML classification.
  //
  //   GEI = log₁₀(1 + (CC × BP × CoverageCv) / (UnitigRatio + ε))
  //
  // The multiplication acts as a soft AND gate: all three numerator terms
  // must be elevated for the product to spike.
  //
  //   Numerator (drivers of chaos): CC × BP × CoverageCv
  //     High CC + high BP = combinatorial path explosion.
  //     CoverageCv distinguishes true polymorphism (uniform coverage, low CV)
  //     from collapsed repeats (reads stacking on hubs, high CV).
  //
  //   Denominator (driver of cleanliness): UnitigRatio + ε
  //     UnitigRatio ≈ 0.95 → divides by ~1 → no amplification.
  //     UnitigRatio ≈ 0.10 → acts as 10× multiplier → penalizes shattered
  //     graphs where nearly every node is a branch point.
  //
  //   Log₁₀ squash:
  //     The raw formula can produce arbitrarily large values. The log₁₀
  //     compresses this to a bounded scale (~0 to 7+).
  //     log₁₀(1 + x) ensures GEI = 0.0 when all inputs are zero.
  //
  // Biological interpretation:
  //   GEI < 0.5:  Pristine — true variant in unique sequence.
  //   GEI 0.5–2:  STR/microsatellite — local loops, anchored flanks.
  //   GEI 2–3:    Complex region — multiple overlapping events.
  //   GEI > 3:    Paralogous collapse — segmental dup / centromeric satellite.
  //
  [[nodiscard]] auto GraphEntanglementIndex() const -> f64 {
    static constexpr f64 EPSILON = 1e-6;
    auto const ccomp = static_cast<f64>(mCyclomaticComplexity);
    auto const bpairs = static_cast<f64>(mNumBranchPoints);
    f64 const raw_chaos = (ccomp * bpairs * mCoverageCv) / (mUnitigRatio + EPSILON);
    return std::log10(1.0 + raw_chaos);
  }

 private:
  usize mCyclomaticComplexity = 0;  ///< M = E - V + 1 (single component)
  usize mNumBranchPoints = 0;       ///< nodes with >= 2 edges in some direction
  usize mMaxSingleDirDegree = 0;    ///< max outgoing edges in any single sign direction
  f64 mUnitigRatio = 0.0;           ///< fraction of nodes with 1-in/1-out
  f64 mCoverageCv = 0.0;            ///< CV of total read support (σ/μ)
  f64 mTipToPathCovRatio = 0.0;     ///< mean tip coverage / mean unitig coverage

  // ComputeGraphComplexity is a friend to populate private fields during
  // O(V+E) graph traversal. All other code uses the public accessors.
  friend auto ComputeGraphComplexity(Graph const& graph, usize component_id) -> GraphComplexity;
};

/// Compute GraphComplexity metrics for a single component.
/// Iterates the Graph's node table filtered by component_id.
/// All metrics are O(V+E) to compute.
[[nodiscard]] auto ComputeGraphComplexity(Graph const& graph, usize component_id)
    -> GraphComplexity;

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_GRAPH_COMPLEXITY_H_
