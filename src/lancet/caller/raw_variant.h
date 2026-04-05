#ifndef SRC_LANCET_CALLER_RAW_VARIANT_H_
#define SRC_LANCET_CALLER_RAW_VARIANT_H_

#include <array>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

#include "absl/container/flat_hash_map.h"
#include "absl/strings/str_cat.h"
#include "lancet/base/sequence_complexity.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/graph_complexity.h"

namespace lancet::caller {

// ============================================================================
// LocusKey: groups RawVariants at the same genomic position for multi-allelic
// VCF output. Two variants with the same LocusKey (chrom, position, ref allele)
// become a single VCF record with comma-separated ALTs.
//
// Example:
//   RawVariant{chr1, 100, A→T}  ─┐
//                                 ├─ same LocusKey → VCF: chr1 100 . A T,G ...
//   RawVariant{chr1, 100, A→G}  ─┘
// ============================================================================
struct LocusKey {
  usize chrom_index = 0;
  usize genome_start1 = 0;
  std::string_view ref_allele;

  friend auto operator==(const LocusKey& lhs, const LocusKey& rhs) -> bool {
    return std::tie(lhs.chrom_index, lhs.genome_start1, lhs.ref_allele) ==
           std::tie(rhs.chrom_index, rhs.genome_start1, rhs.ref_allele);
  }

  friend auto operator<(const LocusKey& lhs, const LocusKey& rhs) -> bool {
    return std::tie(lhs.chrom_index, lhs.genome_start1, lhs.ref_allele) <
           std::tie(rhs.chrom_index, rhs.genome_start1, rhs.ref_allele);
  }

  template <typename H>
  friend auto AbslHashValue(H hash_state, const LocusKey& key) -> H {
    return H::combine(std::move(hash_state), key.chrom_index, key.genome_start1, key.ref_allele);
  }
};

class RawVariant {
 public:
  RawVariant() = default;

  enum class Type : i8 { REF = -1, SNV = 0, INS = 1, DEL = 2, MNP = 3 };
  enum class State : i8 { NONE = -1, SHARED = 0, NORMAL = 1, TUMOR = 2, UNKNOWN = 3 };

  // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
  usize mChromIndex = -1;
  usize mGenomeStart1 = -1;
  i64 mAlleleLength = -1;
  std::string mChromName;
  std::string mRefAllele;
  std::string mAltAllele;

  // haplotype index identifier -> start index of variant in haplotype
  absl::flat_hash_map<usize, usize> mHapStart0Idxs;

  Type mType = Type::REF;

  // ── Graph complexity metrics (ML-ready orthogonal features) ────────────
  // Populated when --enable-graph-complexity-features is set.
  // 3 fields matching the GRAPH_CX VCF INFO tag.
  //
  // Raw topology metrics (CC, BP, EdgeDensity, UnitigRatio, CoverageCv) are
  // mathematically collinear and compressed into the Graph Entanglement Index.
  // Color-based metrics (UnsharedColorRatio, ColorDiscordantBranches) are not
  // topological and are captured by other biologically relevant annotations.
  struct GraphMetrics {
    f64 graph_entanglement_index = 0.0;  ///< GEI: log₁₀(1 + CC×BP×CovCV / UnitigRatio)
    f64 tip_to_path_cov_ratio = 0.0;    ///< assembly tearing: tip cov / unitig cov
    usize max_single_dir_degree = 0;    ///< hub k-mer detection: max outgoing edges

    /// Format as 3 comma-separated values for VCF GRAPH_CX INFO tag.
    [[nodiscard]] auto FormatVcfValue() const -> std::string {
      return absl::StrCat(
          base::FormatComplexityScore(graph_entanglement_index), ",",
          base::FormatComplexityScore(tip_to_path_cov_ratio), ",",
          max_single_dir_degree);
    }
  };

  // ── Annotation fields (mutable — populated post-construction) ─────────
  // These do not participate in btree_set ordering and are annotated after
  // variant discovery, so they are mutable to allow modification through
  // const btree_set iterators without const_cast.

  mutable GraphMetrics mGraphMetrics;

  // ── Sequence complexity (11 ML-ready features) ────────────────────────
  // Populated when --enable-sequence-complexity-features is set.
  // Distilled from raw multi-scale metrics (HRun, Entropy, LongdustQ, TR motifs)
  // into 3 groups: Context (REF), Delta (ALT−REF), TR Motif (ALT).
  mutable base::SequenceComplexity mSeqCx;
  // NOLINTEND(misc-non-private-member-variables-in-classes)

  // Create a key for grouping variants at the same locus into multi-allelic records.
  [[nodiscard]] auto MakeLocusKey() const -> LocusKey {
    return {mChromIndex, mGenomeStart1, mRefAllele};
  }

  template <typename HashState>
  friend auto AbslHashValue(HashState hash_state, const RawVariant& var) -> HashState {
    return HashState::combine(std::move(hash_state), var.mChromIndex, var.mGenomeStart1, var.mAlleleLength,
                              static_cast<i8>(var.mType), var.mChromName, var.mRefAllele, var.mAltAllele);
  }

  friend auto operator==(const RawVariant& lhs, const RawVariant& rhs) -> bool {
    return lhs.mChromIndex == rhs.mChromIndex && lhs.mGenomeStart1 == rhs.mGenomeStart1 &&
           lhs.mAlleleLength == rhs.mAlleleLength && lhs.mType == rhs.mType && lhs.mChromName == rhs.mChromName &&
           lhs.mRefAllele == rhs.mRefAllele && lhs.mAltAllele == rhs.mAltAllele;
  }

  friend auto operator<(const RawVariant& lhs, const RawVariant& rhs) -> bool {
    // NOLINTBEGIN(readability-braces-around-statements)
    if (lhs.mChromIndex != rhs.mChromIndex) return lhs.mChromIndex < rhs.mChromIndex;
    if (lhs.mGenomeStart1 != rhs.mGenomeStart1) return lhs.mGenomeStart1 < rhs.mGenomeStart1;
    if (lhs.mRefAllele != rhs.mRefAllele) return lhs.mRefAllele < rhs.mRefAllele;
    return lhs.mAltAllele < rhs.mAltAllele;
    // NOLINTEND(readability-braces-around-statements)
  }

  friend auto operator!=(const RawVariant& lhs, const RawVariant& rhs) -> bool { return !(rhs == lhs); }
  friend auto operator>(const RawVariant& lhs, const RawVariant& rhs) -> bool { return rhs < lhs; }
  friend auto operator<=(const RawVariant& lhs, const RawVariant& rhs) -> bool { return !(rhs < lhs); }
  friend auto operator>=(const RawVariant& lhs, const RawVariant& rhs) -> bool { return !(lhs < rhs); }
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_RAW_VARIANT_H_
