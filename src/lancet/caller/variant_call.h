#ifndef SRC_LANCET_CALLER_VARIANT_CALL_H_
#define SRC_LANCET_CALLER_VARIANT_CALL_H_

#include <array>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_support.h"
#include "lancet/core/sample_info.h"

namespace lancet::caller {

using VariantID = u64;

// ============================================================================
// VariantCall: a finalized VCF record, potentially multi-allelic.
//
// Multi-allelic support:
//   The constructor accepts a group of RawVariants at the same locus
//   (same chrom, position, ref allele) and merges their per-sample evidence
//   into a unified multi-allelic representation:
//
//     Input:  {chr1:100 A→T, chr1:100 A→G}  (two RawVariants, same locus)
//     Output: VCF record: chr1 100 . A T,G ... (comma-separated ALTs)
//
// FORMAT fields (per-sample):
//   GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:PL:GQ
//
//   GT   - Genotype (e.g., 0/1, 1/2 for multi-allelic)
//   AD   - Number=R: read depth per allele (REF, ALT1, ALT2, ...)
//   ADF  - Number=R: forward strand depth per allele
//   ADR  - Number=R: reverse strand depth per allele
//   DP   - Total read depth
//   RMQ  - Number=R: RMS mapping quality per allele
//   NPBQ - Number=R: normalized posterior base quality per allele (PBQ/N)
//   SB   - Number=1: Strand bias log odds ratio (Haldane-corrected)
//   SCA  - Number=1: Soft Clip Asymmetry (ALT - REF soft-clip fraction)
//   FLD  - Number=1: Fragment Length Delta (|mean ALT isize - mean REF isize|)
//   RPCD - Number=1: Read Position Cohen's D (folded position effect size)
//   BQCD - Number=1: Base Quality Cohen's D (base quality effect size)
//   MQCD - Number=1: Mapping Quality Cohen's D (MAPQ effect size)
//   ASMD - Number=1: Allele-Specific Mismatch Delta (mean ALT NM - mean REF NM)
//   SDFC - Number=1: Site Depth Fold Change (DP / window mean coverage)
//   PRAD - Number=1: Polar Radius log10(1 + sqrt(AD_Ref² + AD_Alt²))
//   PANG - Number=1: Polar Angle atan2(AD_Alt, AD_Ref) in radians
//   PL   - Number=G: Phred-scaled genotype likelihoods
//   GQ   - Genotype quality (second-lowest PL, capped at 99)
// ============================================================================
class VariantCall {
 public:
  using Samples = absl::Span<const core::SampleInfo>;

  /// Feature gate flags — replaces the old string-based AnnotationFeatures.
  struct FeatureFlags {
    bool enable_graph_complexity = false;
    bool enable_sequence_complexity = false;
  };

  // Single-variant (bi-allelic) constructor
  using Supports = absl::flat_hash_map<std::string_view, std::unique_ptr<VariantSupport>>;
  VariantCall(const RawVariant* var, Supports&& supports, Samples samps, FeatureFlags features, f64 window_cov);

  // Multi-allelic constructor: group of variants at the same locus
  using VariantGroup = absl::Span<const RawVariant* const>;
  using SupportsByVariant = absl::flat_hash_map<const RawVariant*, Supports>;
  VariantCall(VariantGroup variants, SupportsByVariant&& all_supports, Samples samps, FeatureFlags features,
              f64 window_cov);

  [[nodiscard]] auto ChromIndex() const -> usize { return mChromIndex; }
  [[nodiscard]] auto ChromName() const -> std::string_view { return mChromName; }
  [[nodiscard]] auto StartPos1() const -> usize { return mStartPos1; }
  [[nodiscard]] auto RefAllele() const -> std::string_view { return mRefAllele; }
  [[nodiscard]] auto AltAlleles() const -> absl::Span<const std::string> { return mAltAlleles; }
  [[nodiscard]] auto NumAltAlleles() const -> usize { return mAltAlleles.size(); }
  [[nodiscard]] auto Length() const -> i64 { return mVariantLength; }
  [[nodiscard]] auto Quality() const -> f64 { return mSiteQuality; }
  [[nodiscard]] auto State() const -> RawVariant::State { return mState; }
  [[nodiscard]] auto Category() const -> RawVariant::Type { return mCategory; }

  [[nodiscard]] auto NumSamples() const -> usize { return mFormatFields.empty() ? 0 : mFormatFields.size() - 1; }
  [[nodiscard]] auto Identifier() const -> VariantID { return mVariantId; }
  [[nodiscard]] auto TotalCoverage() const -> usize { return mTotalSampleCov; }

  /// Returns true if any sample has ALT allele support.
  /// Used by VariantStore to filter out zero-evidence calls.
  [[nodiscard]] auto HasAltSupport() const -> bool { return mHasAltSupport; }

  [[nodiscard]] auto AsVcfRecord() const -> std::string;

  friend auto operator==(const VariantCall& lhs, const VariantCall& rhs) -> bool {
    return lhs.mVariantId == rhs.mVariantId;
  }
  friend auto operator<(const VariantCall& lhs, const VariantCall& rhs) -> bool {
    // NOLINTBEGIN(readability-braces-around-statements)
    if (lhs.mChromIndex != rhs.mChromIndex) return lhs.mChromIndex < rhs.mChromIndex;
    if (lhs.mStartPos1 != rhs.mStartPos1) return lhs.mStartPos1 < rhs.mStartPos1;
    if (lhs.mRefAllele != rhs.mRefAllele) return lhs.mRefAllele < rhs.mRefAllele;
    return lhs.mVariantId < rhs.mVariantId;
    // NOLINTEND(readability-braces-around-statements)
  }

 private:
  u64 mVariantId;
  usize mChromIndex;
  usize mStartPos1;
  usize mTotalSampleCov;
  std::string mChromName;
  std::string mRefAllele;
  std::vector<std::string> mAltAlleles;  // Multi-allelic: [ALT1, ALT2, ...]

  i64 mVariantLength;
  f64 mSiteQuality;
  RawVariant::State mState = RawVariant::State::NONE;
  RawVariant::Type mCategory = RawVariant::Type::REF;
  bool mHasAltSupport = false;  ///< true if any sample has ALT coverage

  std::string mInfoField;

  // ── Graph complexity metrics (from RawVariant, transcribed) ────────────
  RawVariant::GraphMetrics mGraphMetrics;

  // ── Sequence complexity (11 ML-ready features, from RawVariant) ─────────
  base::SequenceComplexity mSeqCx;

  // ── Feature gate flags (from constructor, preserved for BuildInfoField) ─
  FeatureFlags mFeatureFlags;

  /// Window-level BAM coverage used as the SDFC denominator.
  /// Set during construction from ProcessWindow context.
  f64 mWindowCov = 0.0;

  /// Site Depth Fold Change: DP / window mean coverage.
  /// Spikes indicate collapsed paralogous mappings; dips indicate mapping holes.
  [[nodiscard]] auto SiteDepthFoldChange() const -> f64 {
    return mWindowCov > 0.0 ? static_cast<f64>(mTotalSampleCov) / mWindowCov : 1.0;
  }

  std::vector<std::string> mFormatFields;

  // ── Evidence collection (shared by both constructors) ──────────────────
  using PerSampleEvidence = absl::flat_hash_map<const core::SampleInfo, std::unique_ptr<VariantSupport>,
                                                core::SampleInfo::Hash, core::SampleInfo::Equal>;

  /// Common finalization after evidence is assembled: builds FORMAT, state, and INFO fields.
  void Finalize(const PerSampleEvidence& evidence, Samples samps, FeatureFlags features);

  // ── Modular field builders ─────────────────────────────────────────────

  /// Build per-sample FORMAT strings (GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:PL:GQ).
  /// Also computes site quality and total coverage. Returns {alt_in_normal, alt_in_tumor}.
  struct AltPresence {
    bool in_normal = false;
    bool in_tumor = false;
  };

  auto BuildFormatFields(const PerSampleEvidence& evidence, Samples samps, bool tumor_normal_mode) -> AltPresence;

  /// Compute SHARED/NORMAL/TUMOR/UNKNOWN state from ALT presence flags.
  /// In non-tumor-normal mode (i.e. normal-only), state is always UNKNOWN.
  void ComputeState(AltPresence alt_presence, bool tumor_normal_mode);

  /// Assemble the INFO field string (TYPE, LENGTH, optional state prefix, optional complexity annotations).
  void BuildInfoField(bool tumor_normal_mode, FeatureFlags features);

  // ── Static helpers ─────────────────────────────────────────────────────

  /// Convert a GL index back to genotype string (e.g., GL=4 with k=3 → "1/2")
  [[nodiscard]] static auto GenotypeFromGLIndex(usize gl_index, usize num_alleles) -> std::string;

  /// Somatic log odds ratio: tumor ALT enrichment vs normal.
  [[nodiscard]] static auto SomaticLogOddsRatio(const core::SampleInfo& curr,
                                                const PerSampleEvidence& supports) -> f64;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_CALL_H_
