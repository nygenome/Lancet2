#include "lancet/core/variant_annotator.h"

#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_set.h"
#include "lancet/cbdg/graph_complexity.h"

#include "absl/types/span.h"

#include <algorithm>
#include <limits>
#include <string>

namespace lancet::core {

// ============================================================================
// Constructor
// ============================================================================

VariantAnnotator::VariantAnnotator(f64 gc_frac) : mSeqCxScorer(gc_frac) {}

// ============================================================================
// AnnotateSequenceComplexity: produce 11-feature SequenceComplexity per variant
//
// For each variant, the scorer needs both REF and ALT haplotypes:
//   - REF = haplotypes[0] (reference haplotype)
//   - ALT = each haplotype carrying the variant (max across haplotypes)
//
// The scorer internally computes context (from REF), deltas (ALT−REF),
// and TR motif features (from ALT), then distills into the 11-feature output.
//
// Multi-haplotype merging: if a variant appears on multiple ALT haplotypes,
// we score each one and MergeMax (element-wise pessimistic worst-case).
//
// If a variant only exists on the REF haplotype (no ALT haps), we score
// REF vs REF to populate context features (deltas will be zero).
// ============================================================================
// NOLINTNEXTLINE(readability-function-size)
void VariantAnnotator::AnnotateSequenceComplexity(caller::VariantSet const& vset,
                                                  absl::Span<std::string const> haplotypes) const {
  static constexpr usize REF_HAP_IDX = 0;

  for (auto const& var : vset) {
    auto const ref_pos = var.mLocalRefStart0Idx;
    if (ref_pos == std::numeric_limits<usize>::max() || REF_HAP_IDX >= haplotypes.size()) {
      continue;
    }

    auto const& ref_hap = haplotypes[REF_HAP_IDX];
    auto const ref_len = var.mRefAllele.size();

    bool scored_any_alt = false;

    // Score each haplotype backing every distinct ALT, dynamically merging the worst-case
    // complexity
    for (auto const& alt : var.mAlts) {
      auto const alt_len = std::max(ref_len, alt.mSequence.size());

      for (auto const& [hap_idx, hap_pos] : alt.mLocalHapStart0Idxs) {
        if (hap_idx >= haplotypes.size() || hap_idx == REF_HAP_IDX) continue;

        auto cplx =
            mSeqCxScorer.Score(ref_hap, ref_pos, ref_len, haplotypes[hap_idx], hap_pos, alt_len);
        var.mSeqCx.MergeMax(cplx);
        scored_any_alt = true;
      }
    }

    // If no ALT haplotypes found, score REF vs REF to populate context features
    if (!scored_any_alt) {
      var.mSeqCx = mSeqCxScorer.Score(ref_hap, ref_pos, ref_len, ref_hap, ref_pos, ref_len);
    }
  }
}

// ============================================================================
// AnnotateGraphComplexity: populate graph complexity metrics on variants
//
// All variants from the same component share the same graph complexity.
// This is by design: graph complexity is a property of the assembly
// component, not individual variants.
//
// Raw topology metrics (CC, BP, EdgeDensity, UnitigRatio, CoverageCv) are
// compressed into the Graph Entanglement Index (GEI) to eliminate collinearity.
// Only orthogonal features are preserved: GEI, TipToPathCovRatio, MaxDegree.
// ============================================================================
void VariantAnnotator::AnnotateGraphComplexity(caller::VariantSet const& vset,
                                               cbdg::GraphComplexity const& component_cx) {
  auto const metrics = caller::RawVariant::GraphMetrics{
      .mGraphEntanglementIndex = component_cx.GraphEntanglementIndex(),
      .mTipToPathCovRatio = component_cx.TipToPathCovRatio(),
      .mMaxSingleDirDegree = component_cx.MaxSingleDirDegree(),
  };

  std::ranges::for_each(vset, [&metrics](auto const& var) -> void { var.mGraphMetrics = metrics; });
}

}  // namespace lancet::core
