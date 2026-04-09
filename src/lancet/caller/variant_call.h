#ifndef SRC_LANCET_CALLER_VARIANT_CALL_H_
#define SRC_LANCET_CALLER_VARIANT_CALL_H_

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
//   The constructor natively ingests a single fully multiallelic RawVariant
//   directly emitted from the Genotyper module. It unpacks the pre-mapped
//   SupportArray layout directly into a unified multidimensional VCF trace:
//
//     Input:  {chr1:100 A→T, A→G}  (One Natively Multiallelic RawVariant)
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

  // Native multi-allelic constructor
  using SupportsByVariant = absl::flat_hash_map<const RawVariant*, SupportArray>;
  VariantCall(const RawVariant* var, const SupportsByVariant& all_supports, Samples samps, FeatureFlags features, f64 window_cov);

  [[nodiscard]] auto ChromIndex() const -> usize { return mChromIndex; }
  [[nodiscard]] auto ChromName() const -> std::string_view { return mChromName; }
  [[nodiscard]] auto StartPos1() const -> usize { return mStartPos1; }
  [[nodiscard]] auto RefAllele() const -> std::string_view { return mRefAllele; }
  [[nodiscard]] auto AltAlleles() const -> absl::Span<const std::string> { return mAltAlleles; }
  [[nodiscard]] auto NumAltAlleles() const -> usize { return mAltAlleles.size(); }
  [[nodiscard]] auto VariantLengths() const -> absl::Span<const i64> { return mVariantLengths; }
  [[nodiscard]] auto Quality() const -> f64 { return mSiteQuality; }
  [[nodiscard]] auto State() const -> RawVariant::State { return mState; }
  [[nodiscard]] auto Categories() const -> absl::Span<const RawVariant::Type> { return mCategories; }
  [[nodiscard]] auto IsMultiallelic() const -> bool { return mIsMultiallelic; }

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
    return lhs.mTotalSampleCov > rhs.mTotalSampleCov;
    // NOLINTEND(readability-braces-around-statements)
  }

 private:
  // ── 8B Alignment ────────────────────────────────────────────────────────
  u64 mVariantId;
  usize mChromIndex;
  usize mStartPos1;
  usize mTotalSampleCov;
  f64 mSiteQuality;
  f64 mWindowCov;

  std::string mChromName;
  std::string mRefAllele;
  std::string mInfoField;
  
  std::vector<i64> mVariantLengths;
  std::vector<std::string> mAltAlleles;
  std::vector<RawVariant::Type> mCategories;
  std::vector<std::string> mFormatFields;

  // ── Sequence complexity (11 ML-ready features, from RawVariant) ─────────
  // ── Graph complexity metrics (from RawVariant, transcribed) ────────────
  base::SequenceComplexity mSeqCx;
  RawVariant::GraphMetrics mGraphCx;

  // ── 4B/2B Alignment ─────────────────────────────────────────────────────
  RawVariant::State mState = RawVariant::State::NONE;

  // ── 1B Alignment ────────────────────────────────────────────────────────
  FeatureFlags mFeatureFlags;
  bool mIsMultiallelic = false;
  bool mHasAltSupport = false;

  /// Site Depth Fold Change: DP / window mean coverage.
  /// Spikes indicate collapsed paralogous mappings; dips indicate mapping holes.
  [[nodiscard]] auto SiteDepthFoldChange() const -> f64 {
    return mWindowCov > 0.0 ? static_cast<f64>(mTotalSampleCov) / mWindowCov : 1.0;
  }

  // ── Evidence collection (shared by both constructors) ──────────────────

  /// Common finalization after evidence is assembled: builds FORMAT, state, and INFO fields.
  void Finalize(const SupportArray& evidence, Samples samps, FeatureFlags features);

  // ── Modular field builders ─────────────────────────────────────────────

  struct AltPresence {
    bool in_normal = false;
    bool in_tumor = false;
  };
  
  struct PerAlleleMetrics {
    std::string ad;
    std::string adf;
    std::string adr;
    std::string rmq;
    std::string npbq;
  };

  /// Build per-sample FORMAT strings (GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:PL:GQ).
  /// Also computes site quality and total coverage. Returns {alt_in_normal, alt_in_tumor}.
  auto BuildFormatFields(const SupportArray& evidence, Samples samps, bool tumor_normal_mode) -> AltPresence;

  [[nodiscard]] auto BuildGenotype(absl::Span<const int> pls, usize num_alleles) const -> std::string;
  
  /// Convert a GL index back to genotype string (e.g., GL=4 with k=3 → "1/2")
  [[nodiscard]] static auto GenotypeFromGLIndex(usize gl_index, usize num_alleles) -> std::string;

  void UpdateSiteQuality(const core::SampleInfo& sinfo, const VariantSupport* support, const SupportArray& evidence, Samples samps, bool tumor_normal_mode, absl::Span<const int> pls);
  
  /// Somatic log odds ratio: tumor ALT enrichment vs normal.
  [[nodiscard]] static auto SomaticLogOddsRatio(const core::SampleInfo& curr, const SupportArray& supports, Samples samps) -> f64;

  void TrackAltPresence(const VariantSupport* support, const core::SampleInfo& sinfo, bool tumor_normal_mode, AltPresence& alt_presence);
  
  /// Generate the Number=R scalar metrics natively iterating across each variant allele index.
  [[nodiscard]] auto BuildPerAlleleMetrics(const VariantSupport* support, usize num_alleles) const -> PerAlleleMetrics;

  /// Compute SHARED/NORMAL/TUMOR/UNKNOWN state from ALT presence flags.
  /// In non-tumor-normal mode (i.e. normal-only), state is always UNKNOWN.
  void ComputeState(AltPresence alt_presence, bool tumor_normal_mode);

  /// Assemble the INFO field string (TYPE, LENGTH, optional state prefix, optional complexity annotations).
  void BuildInfoField(bool tumor_normal_mode, FeatureFlags features);
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_CALL_H_
