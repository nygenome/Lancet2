#ifndef SRC_LANCET_CORE_VARIANT_ANNOTATOR_H_
#define SRC_LANCET_CORE_VARIANT_ANNOTATOR_H_

#include <string>

#include "absl/types/span.h"
#include "lancet/base/sequence_complexity.h"
#include "lancet/base/types.h"
#include "lancet/caller/variant_set.h"
#include "lancet/cbdg/graph_complexity.h"

namespace lancet::core {

// ============================================================================
// VariantAnnotator — annotates RawVariants with complexity features
//
// Separated from VariantBuilder to isolate annotation logic (which produces
// ML-ready features) from variant discovery logic (which walks the graph).
//
// Owns the SequenceComplexityScorer and provides methods to annotate
// variant sets with sequence complexity and graph complexity metrics.
// ============================================================================
class VariantAnnotator {
 public:
  /// @param gc_frac  Global background GC fraction for LongdustQ scoring.
  ///                 Default: 0.41 (human genome-wide average).
  explicit VariantAnnotator(f64 gc_frac = 0.41);

  /// Annotate all variants in `vset` with 11-feature sequence complexity.
  /// Scores each variant against REF and ALT haplotypes, merging across
  /// multiple ALT haplotypes via element-wise max (pessimistic worst-case).
  void AnnotateSequenceComplexity(const caller::VariantSet& vset,
                                  absl::Span<const std::string> haplotypes) const;

  /// Annotate all variants in `vset` with graph complexity metrics.
  /// All variants from the same component share the same graph complexity
  /// (graph complexity is a property of the assembly component, not
  /// individual variants).
  static void AnnotateGraphComplexity(const caller::VariantSet& vset,
                                      const cbdg::GraphComplexity& component_cx);

 private:
  /// Sequence complexity scorer — produces 11-feature SequenceComplexity per variant.
  /// Owns its LongdustQScorer instances. Constructed once with gc_frac.
  base::SequenceComplexityScorer mSeqCxScorer;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_VARIANT_ANNOTATOR_H_
