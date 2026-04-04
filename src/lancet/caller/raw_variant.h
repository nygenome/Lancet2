#ifndef SRC_LANCET_CALLER_RAW_VARIANT_H_
#define SRC_LANCET_CALLER_RAW_VARIANT_H_

#include <array>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

#include "absl/container/flat_hash_map.h"
#include "lancet/base/lcr_scorer.h"
#include "lancet/base/types.h"

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
  enum class State : i8 { NONE = -1, SHARED = 0, NORMAL = 1, TUMOR = 2 };

  // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
  usize mChromIndex = -1;
  usize mGenomeStart1 = -1;
  i64 mAlleleLength = -1;
  Type mType = Type::REF;
  std::string mChromName;
  std::string mRefAllele;
  std::string mAltAllele;

  // haplotype index identifier -> start index of variant in haplotype
  absl::flat_hash_map<usize, usize> mHapStart0Idxs;

  // ALT_LCR: multi-scale low-complexity scores centered on the ALT allele position.
  // Scored across all haplotypes (REF + ALTs) carrying this variant. Max across haplotypes.
  // Indexed by scale: [5bp, 10bp, 50bp, 100bp, full_haplotype].
  // Higher values indicate more repetitive context. Score >= 0.6 ≈ LCR threshold.
  std::array<f64, base::NUM_LCR_SCALES> mAltLcrScores = base::DEFAULT_LCR_SCORES;

  // REF_LCR: multi-scale low-complexity scores centered on the REF allele position.
  // Scored exclusively on the reference haplotype (comp_haps[0]).
  // Indexed by scale: [5bp, 10bp, 50bp, 100bp, full_haplotype].
  std::array<f64, base::NUM_LCR_SCALES> mRefLcrScores = base::DEFAULT_LCR_SCORES;
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
