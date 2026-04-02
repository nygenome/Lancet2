#ifndef SRC_LANCET_CALLER_VARIANT_CALL_H_
#define SRC_LANCET_CALLER_VARIANT_CALL_H_

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
//   GT:AD:ADF:ADR:DP:RMQ:PBQ:SB:PL:GQ
//
//   GT  - Genotype (e.g., 0/1, 1/2 for multi-allelic)
//   AD  - Number=R: read depth per allele (REF, ALT1, ALT2, ...)
//   ADF - Number=R: forward strand depth per allele
//   ADR - Number=R: reverse strand depth per allele
//   DP  - Total read depth
//   RMQ - Number=R: RMS mapping quality per allele
//   PBQ - Number=R: posterior base quality per allele
//   SB  - Number=R: strand bias ratio per allele
//   PL  - Number=G: Phred-scaled genotype likelihoods
//   GQ  - Genotype quality (second-lowest PL, capped at 99)
// ============================================================================
class VariantCall {
 public:
  using Samples = absl::Span<const core::SampleInfo>;

  // Single-variant (bi-allelic) constructor
  using Supports = absl::flat_hash_map<std::string_view, std::unique_ptr<VariantSupport>>;
  VariantCall(const RawVariant* var, Supports&& supports, Samples samps);

  // Multi-allelic constructor: group of variants at the same locus
  using VariantGroup = absl::Span<const RawVariant* const>;
  using SupportsByVariant = absl::flat_hash_map<const RawVariant*, Supports>;
  VariantCall(VariantGroup variants, SupportsByVariant&& all_supports, Samples samps);

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
  RawVariant::State mState;
  RawVariant::Type mCategory;

  std::string mInfoField;
  std::vector<std::string> mFormatFields;

  // Convert a GL index back to genotype string (e.g., GL=4 with k=3 → "1/2")
  [[nodiscard]] static auto GenotypeFromGLIndex(usize gl_index, usize num_alleles) -> std::string;

  using PerSampleEvidence = absl::flat_hash_map<const core::SampleInfo, std::unique_ptr<VariantSupport>,
                                                core::SampleInfo::Hash, core::SampleInfo::Equal>;

  void BuildFormatFields(const PerSampleEvidence& per_sample_evidence, Samples samps, bool germline_mode);

  [[nodiscard]] static auto SomaticFisherScore(const core::SampleInfo& curr,
                                                const PerSampleEvidence& supports) -> f64;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_CALL_H_
