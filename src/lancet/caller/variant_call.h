#ifndef SRC_LANCET_CALLER_VARIANT_CALL_H_
#define SRC_LANCET_CALLER_VARIANT_CALL_H_

#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_support.h"
#include "lancet/core/sample_info.h"

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"

#include <string>
#include <string_view>
#include <vector>

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
  using Samples = absl::Span<core::SampleInfo const>;

  /// Feature gate flags — replaces the old string-based AnnotationFeatures.
  struct FeatureFlags {
    bool mEnableGraphComplexity = false;
    bool mEnableSequenceComplexity = false;
  };

  // Native multi-allelic constructor
  using SupportsByVariant = absl::flat_hash_map<RawVariant const*, SupportArray>;
  VariantCall(RawVariant const* var, SupportsByVariant const& all_supports, Samples samps,
              FeatureFlags features, f64 window_cov);

  [[nodiscard]] auto ChromIndex() const -> usize { return mChromIndex; }
  [[nodiscard]] auto ChromName() const -> std::string_view { return mChromName; }
  [[nodiscard]] auto StartPos1() const -> usize { return mStartPos1; }
  [[nodiscard]] auto RefAllele() const -> std::string_view { return mRefAllele; }
  [[nodiscard]] auto AltAlleles() const -> absl::Span<std::string const> { return mAltAlleles; }
  [[nodiscard]] auto NumAltAlleles() const -> usize { return mAltAlleles.size(); }
  [[nodiscard]] auto VariantLengths() const -> absl::Span<i64 const> { return mVariantLengths; }
  [[nodiscard]] auto Quality() const -> f64 { return mSiteQuality; }
  [[nodiscard]] auto State() const -> RawVariant::State { return mState; }
  [[nodiscard]] auto Categories() const -> absl::Span<RawVariant::Type const> {
    return mCategories;
  }
  [[nodiscard]] auto IsMultiallelic() const -> bool { return mIsMultiallelic; }

  [[nodiscard]] auto NumSamples() const -> usize {
    return mFormatFields.empty() ? 0 : mFormatFields.size() - 1;
  }
  [[nodiscard]] auto Identifier() const -> VariantID { return mVariantId; }
  [[nodiscard]] auto TotalCoverage() const -> usize { return mTotalSampleCov; }

  /// Returns true if any sample has ALT allele support.
  /// Used by VariantStore to filter out zero-evidence calls.
  [[nodiscard]] auto HasAltSupport() const -> bool { return mHasAltSupport; }

  [[nodiscard]] auto AsVcfRecord() const -> std::string;

  // ── VARIANT IDENTITY & ORDERING DESIGN ──────────────────────────────────
  // Identity (operator==) and ordering (operator<) both use the same
  // conceptual key: CHROM + POS + REF (locus-level, ALTs excluded).
  //
  // WHY NO ALTs IN IDENTITY:
  //   Overlapping genomic windows independently assemble the same locus.
  //   Window A might produce chr1:100 A→T at 80x, while Window B produces
  //   chr1:100 A→T,G at 120x. These are the same locus — the higher-coverage
  //   window assembled a more complete multi-allelic picture. With ALTs in the
  //   hash, they'd be treated as *different* variants, producing duplicate VCF
  //   records at the same locus (poor practice per VCF v4.5 spec). With
  //   CHROM+POS+REF identity, VariantStore dedup keeps only the higher-coverage
  //   call, ensuring at most one VCF record per locus.
  //
  // DOWNSTREAM INVARIANT (VariantStore):
  //   VariantStore uses Identifier() (== mVariantId) as flat_hash_map key.
  //   When a duplicate is found (same CHROM+POS+REF), the variant with higher
  //   TotalCoverage() replaces the existing one. This means after dedup, at
  //   most one variant per locus exists, so the TotalCov and mVariantId
  //   tiebreakers in operator< are purely for mathematical completeness.
  //
  // WARNING: Do NOT add ALT alleles to HashRawVariant or change the fields
  //   used by operator== / operator< without updating VariantStore's dedup
  //   logic. Mismatched identity semantics will produce duplicate VCF records
  //   and break tabix indexing.
  // ────────────────────────────────────────────────────────────────────────
  friend auto operator==(VariantCall const& lhs, VariantCall const& rhs) -> bool {
    return lhs.mVariantId == rhs.mVariantId;
  }
  friend auto operator<(VariantCall const& lhs, VariantCall const& rhs) -> bool {
    if (lhs.mChromIndex != rhs.mChromIndex) return lhs.mChromIndex < rhs.mChromIndex;
    if (lhs.mStartPos1 != rhs.mStartPos1) return lhs.mStartPos1 < rhs.mStartPos1;
    if (lhs.mRefAllele != rhs.mRefAllele) return lhs.mRefAllele < rhs.mRefAllele;

    if (lhs.mTotalSampleCov != rhs.mTotalSampleCov) {
      return lhs.mTotalSampleCov < rhs.mTotalSampleCov;
    }

    return lhs.mVariantId < rhs.mVariantId;
  }

 private:
  // ── 8B Alignment ────────────────────────────────────────────────────────
  u64 mVariantId;
  usize mChromIndex;
  usize mStartPos1;
  usize mTotalSampleCov{0};
  f64 mSiteQuality{0};
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
  void Finalize(SupportArray const& evidence, Samples samps, FeatureFlags features);

  // ── Modular field builders ─────────────────────────────────────────────

  struct AltPresence {
    bool mInNormal = false;
    bool mInTumor = false;
  };

  struct PerAlleleMetrics {
    std::string mAd;
    std::string mAdf;
    std::string mAdr;
    std::string mRmq;
    std::string mNpbq;
  };

  /// Build per-sample FORMAT strings
  /// (GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:PL:GQ).
  /// Also computes site quality and total coverage. Returns {alt_in_normal, alt_in_tumor}.
  auto BuildFormatFields(SupportArray const& evidence, Samples samps, bool tumor_normal_mode)
      -> AltPresence;

  [[nodiscard]] static auto BuildGenotype(absl::Span<int const> pls, usize num_alleles)
      -> std::string;

  /// Convert a GL index back to genotype string (e.g., GL=4 with k=3 → "1/2")
  [[nodiscard]] static auto GenotypeFromGLIndex(usize gl_index, usize num_alleles) -> std::string;

  void UpdateSiteQuality(core::SampleInfo const& sinfo, VariantSupport const* support,
                         SupportArray const& evidence, Samples samps, bool tumor_normal_mode,
                         absl::Span<int const> pls);

  /// Somatic log odds ratio: tumor ALT enrichment vs normal.
  [[nodiscard]] static auto SomaticLogOddsRatio(core::SampleInfo const& curr,
                                                SupportArray const& supports, Samples samps) -> f64;

  void TrackAltPresence(VariantSupport const* support, core::SampleInfo const& sinfo,
                        bool tumor_normal_mode, AltPresence& alt_presence);

  /// Generate the Number=R scalar metrics natively iterating across each variant allele index.
  [[nodiscard]] static auto BuildPerAlleleMetrics(VariantSupport const* support, usize num_alleles)
      -> PerAlleleMetrics;

  /// Compute SHARED/NORMAL/TUMOR/UNKNOWN state from ALT presence flags.
  /// In non-tumor-normal mode (i.e. normal-only), state is always UNKNOWN.
  void ComputeState(AltPresence alt_presence, bool tumor_normal_mode);

  /// Assemble the INFO field string (TYPE, LENGTH, optional state prefix, optional complexity
  /// annotations).
  void BuildInfoField(bool tumor_normal_mode, FeatureFlags features);
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_CALL_H_
