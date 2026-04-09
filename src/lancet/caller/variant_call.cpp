#include "lancet/caller/variant_call.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/hash/hash.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "lancet/base/assert.h"
#include "lancet/base/longdust_scorer.h"
#include "lancet/base/polar_coords.h"
#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_support.h"
#include "lancet/cbdg/label.h"
#include "spdlog/fmt/bundled/core.h"

namespace {

[[nodiscard]] inline auto HashRawVariant(const lancet::caller::RawVariant* var) -> u64 {
  std::size_t h = absl::HashOf(var->mChromName, var->mGenomeChromPos1, var->mRefAllele);
  std::for_each(var->mAlts.cbegin(), var->mAlts.cend(), [&h](const auto& alt) {
    h = absl::HashOf(h, alt.mSequence, alt.mLength, alt.mType);
  });
  return static_cast<u64>(h);
}

}  // namespace

namespace lancet::caller {

// ============================================================================
// Multi-allelic constructor
//
// Unspools a natively multi-allelic RawVariant exactly into the VCF array mapping!
// ============================================================================
VariantCall::VariantCall(const RawVariant* var, const SupportsByVariant& all_supports, Samples samps,
                         FeatureFlags features, const f64 window_cov)
    : mVariantId(HashRawVariant(var)),
      mChromIndex(var->mChromIndex), mStartPos1(var->mGenomeChromPos1),
      mTotalSampleCov(0), mSiteQuality(0), mWindowCov(window_cov),
      mChromName(var->mChromName), mRefAllele(var->mRefAllele),
      mFeatureFlags(features) {
          
  const auto num_alts = var->mAlts.size();
  mAltAlleles.reserve(num_alts);
  mCategories.reserve(num_alts);
  mVariantLengths.reserve(num_alts);

  std::for_each(var->mAlts.cbegin(), var->mAlts.cend(), [this](const auto& alt) {
    mAltAlleles.push_back(alt.mSequence);
    mCategories.push_back(alt.mType);
    mVariantLengths.push_back(alt.mLength);
  });
  
  mIsMultiallelic = (mAltAlleles.size() > 1);
  mSeqCx = var->mSeqCx;
  mGraphCx = var->mGraphMetrics;

  if (const auto var_it = all_supports.find(var); var_it != all_supports.end()) {
    Finalize(var_it->second, samps, features);
  } else {
    Finalize(SupportArray(), samps, features);
  }
}

// ============================================================================
// Finalize: common post-construction steps for both constructors.
//
// Determines tumor-normal (somatic) vs normal-only mode, then delegates to
// three focused methods:
//   1. BuildFormatFields  — per-sample FORMAT strings + site quality
//   2. ComputeState       — SHARED/NORMAL/TUMOR classification
//   3. BuildInfoField     — INFO string assembly
// ============================================================================
void VariantCall::Finalize(const SupportArray& evidence, Samples samps, FeatureFlags features) {
  static const auto is_tumor = [](const auto& s) -> bool { return s.TagKind() == cbdg::Label::TUMOR; };
  static const auto is_normal = [](const auto& s) -> bool { return s.TagKind() == cbdg::Label::NORMAL; };
  const auto tumor_normal_mode = std::ranges::any_of(samps, is_tumor) && std::ranges::any_of(samps, is_normal);

  const auto alt_presence = BuildFormatFields(evidence, samps, tumor_normal_mode);
  ComputeState(alt_presence, tumor_normal_mode);
  BuildInfoField(tumor_normal_mode, features);
}

// ============================================================================
// BuildFormatFields: per-sample FORMAT strings and site quality.
//
// FORMAT: GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:PL:GQ
//
//   GT   - Genotype derived from minimum PL
//   AD   - Number=R: allele depths (REF, ALT1, ALT2, ...)
//   ADF  - Number=R: forward strand allele depths
//   ADR  - Number=R: reverse strand allele depths
//   DP   - Total depth
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
//   GQ   - Genotype quality (GATK: second-lowest PL, capped at 99)
// ============================================================================
auto VariantCall::BuildFormatFields(const SupportArray& evidence, Samples samps, const bool tumor_normal_mode) -> AltPresence {
  const auto num_alleles = mAltAlleles.size() + 1;  // +1 for REF
  AltPresence alt_presence;

  mFormatFields.reserve(samps.size() + 1);
  mFormatFields.emplace_back("GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:PL:GQ");

  for (const auto& sinfo : samps) {
    const auto* support = evidence.Find(sinfo.SampleName());
    if (support == nullptr) continue;

    mTotalSampleCov += support->TotalSampleCov();
    
    const auto pls = support->ComputePLs();
    const auto gq_value = VariantSupport::ComputeGQ(pls);
    const auto genotype = BuildGenotype(pls, num_alleles);
    const auto allele_metrics = BuildPerAlleleMetrics(support, num_alleles);

    UpdateSiteQuality(sinfo, support, evidence, samps, tumor_normal_mode, pls);
    TrackAltPresence(support, sinfo, tumor_normal_mode, alt_presence);

    // Strand bias log odds ratio (Number=1, per-sample)
    const auto sb_val = fmt::format("{:.3f}", support->StrandBiasLogOR());

    // Alignment-derived per-sample annotations (coverage-normalized effect sizes)
    const auto sca_val = fmt::format("{:.4f}", support->SoftClipAsymmetry());
    const auto fld_val = fmt::format("{:.1f}", support->FragLengthDelta());
    const auto rpcd_val = fmt::format("{:.4f}", support->ReadPosCohenD());
    const auto bqcd_val = fmt::format("{:.4f}", support->BaseQualCohenD());
    const auto mqcd_val = fmt::format("{:.4f}", support->MappingQualCohenD());
    const auto asmd_val = fmt::format("{:.3f}", support->AlleleMismatchDelta());
    const auto sdfc_val = fmt::format("{:.2f}", SiteDepthFoldChange());

    // Polar coordinate features for ML variant classification
    // PRAD/PANG orthogonalize allele identity from depth (see polar_coords.h)
    const auto ad_ref = static_cast<f64>(support->TotalRefCov());
    const auto ad_alt = static_cast<f64>(support->TotalAltCov());
    const auto prad_val = fmt::format("{:.4f}", base::PolarRadius(ad_ref, ad_alt));
    const auto pang_val = fmt::format("{:.4f}", base::PolarAngle(ad_alt, ad_ref));

    const std::string pl_val = pls.empty() ? "." : absl::StrJoin(pls, ",");

    mFormatFields.emplace_back(fmt::format(
        "{GT}:{AD}:{ADF}:{ADR}:{DP}:{RMQ}:{NPBQ}:{SB}:{SCA}:{FLD}:{RPCD}:{BQCD}:{MQCD}:{ASMD}:{SDFC}:{PRAD}:{PANG}:{PL}:{GQ}",
        fmt::arg("GT", genotype), fmt::arg("AD", allele_metrics.ad), fmt::arg("ADF", allele_metrics.adf), fmt::arg("ADR", allele_metrics.adr),
        fmt::arg("DP", support->TotalSampleCov()), fmt::arg("RMQ", allele_metrics.rmq), fmt::arg("NPBQ", allele_metrics.npbq),
        fmt::arg("SB", sb_val), fmt::arg("SCA", sca_val), fmt::arg("FLD", fld_val), fmt::arg("RPCD", rpcd_val),
        fmt::arg("BQCD", bqcd_val), fmt::arg("MQCD", mqcd_val), fmt::arg("ASMD", asmd_val), fmt::arg("SDFC", sdfc_val),
        fmt::arg("PRAD", prad_val), fmt::arg("PANG", pang_val), fmt::arg("PL", pl_val), fmt::arg("GQ", gq_value))
      );
  }

  return alt_presence;
}

auto VariantCall::BuildGenotype(absl::Span<const int> pls, usize num_alleles) const -> std::string {
  usize best_gt_idx = 0;
  if (!pls.empty()) {
    best_gt_idx = static_cast<usize>(std::distance(pls.cbegin(), std::min_element(pls.cbegin(), pls.cend())));
  }
  return GenotypeFromGLIndex(best_gt_idx, num_alleles);
}

// ============================================================================
// GenotypeFromGLIndex: convert a GL index back to a genotype string.
//
// VCF 4.3 §1.6.2 defines the GL index for a diploid genotype (i, j) as:
//
//   GL_index = j * (j + 1) / 2 + i     where i ≤ j
//
// This maps genotypes to a flat array indexed by triangular numbers:
//
//   Genotype:  0/0  0/1  1/1  0/2  1/2  2/2  0/3  1/3  2/3  3/3
//   GL index:   0    1    2    3    4    5    6    7    8    9
//                    └ j=1 ┘    └── j=2 ──┘    └─── j=3 ───┘
//
// Each "row" j starts at T(j) = j*(j+1)/2 and contains i = 0..j.
//
// To invert, we find which row j the GL index falls in, then recover i.
// We use the same integer-only algorithm as htslib's bcf_gt2alleles()
// (htslib/vcf.h): walk triangular numbers T(1), T(2), ... until T(dk) ≥ igt.
//
//   k  = running triangular number T(dk-1)
//   dk = next row size (starts at 1, incremented each step)
//
//   After the loop: dk-1 = j, and i = igt - T(j) = igt - k + j
//
// Example: gl_index = 4
//   dk=1, k=0 → k<4, dk=2, k=2 → k<4, dk=3, k=5 → k≥4, stop
//   j = dk-1 = 2,  i = 4 - 5 + 2 = 1  → "1/2" ✓
// ============================================================================
auto VariantCall::GenotypeFromGLIndex(const usize gl_index, const usize num_alleles) -> std::string {
  // Integer-only triangular number walk (matches htslib bcf_gt2alleles)
  usize k = 0;
  usize dk = 1;
  while (k < gl_index) {
    dk++;
    k += dk;
  }
  const auto j = dk - 1;
  const auto i = gl_index - k + j;

  // This should never trigger: gl_index comes from std::min_element over the
  // PL vector (size = num_alleles*(num_alleles+1)/2), so the inverted (i,j)
  // will always satisfy i <= j < num_alleles. If it does fire, ComputePLs or
  // its caller has a bug — don't silently return "0/0" and mask it.
  LANCET_ASSERT(i < num_alleles && j < num_alleles);
  return fmt::format("{}/{}", i, j);
}

void VariantCall::UpdateSiteQuality(const core::SampleInfo& sinfo, const VariantSupport* support,
                                    const SupportArray& evidence, Samples samps,
                                    bool tumor_normal_mode, absl::Span<const int> pls) {
  const auto somatic_lor = tumor_normal_mode ? SomaticLogOddsRatio(sinfo, evidence, samps) : 0.0;
  const auto ref_hom_pl = pls.empty() ? 0 : pls[0];
  const auto sample_dp = support->TotalSampleCov();
  const auto per_read_qual = (!pls.empty() && sample_dp > 0) ? std::min(static_cast<f64>(ref_hom_pl) / static_cast<f64>(sample_dp), 10.0) : 0.0;
  mSiteQuality = std::max(mSiteQuality, tumor_normal_mode ? somatic_lor : per_read_qual);
}

// ============================================================================
// SomaticLogOddsRatio: somatic variant evidence as a coverage-invariant
// log odds ratio. Compares ALT/REF allele counts between the current tumor
// sample and the average across normal samples.
//
//   SOLOR = ln( ((tmr_alt+1)(nml_ref+1)) / ((tmr_ref+1)(nml_alt+1)) )
//
// Haldane correction (+1) handles zero-count edge cases without sentinels.
// A clean somatic produces SOLOR ≈ 5; germline produces SOLOR ≈ 0.
// Coverage-stable: once normal VAF stabilizes (≥60×), SOLOR varies < 3%.
// ============================================================================
auto VariantCall::SomaticLogOddsRatio(const core::SampleInfo& curr, const SupportArray& supports, Samples samps) -> f64 {
  if (curr.TagKind() != cbdg::Label::TUMOR) {
    return 0.0;
  }

  const auto* tmr_evidence = supports.Find(curr.SampleName());
  // Haldane correction (+1 to all cells) avoids mathematical undefined zero-division
  const f64 ta = (tmr_evidence != nullptr ? static_cast<f64>(tmr_evidence->TotalAltCov()) : 0.0) + 1.0;
  const f64 tr = (tmr_evidence != nullptr ? static_cast<f64>(tmr_evidence->TotalRefCov()) : 0.0) + 1.0;

  u32 sum_na = 0;
  u32 sum_nr = 0;
  u32 count_nml = 0;

  for (const auto& sinfo : samps) {
    if (sinfo.TagKind() == cbdg::Label::NORMAL) {
      if (const auto* evidence = supports.Find(sinfo.SampleName()); evidence != nullptr) {
        sum_na += evidence->TotalAltCov();
        sum_nr += evidence->TotalRefCov();
        count_nml++;
      }
    }
  }

  const f64 mean_na = count_nml > 0 ? static_cast<f64>(sum_na) / count_nml : 0.0;
  const f64 mean_nr = count_nml > 0 ? static_cast<f64>(sum_nr) / count_nml : 0.0;

  const f64 na = mean_na + 1.0;
  const f64 nr = mean_nr + 1.0;

  return std::log((ta * nr) / (tr * na));
}

void VariantCall::TrackAltPresence(const VariantSupport* support, const core::SampleInfo& sinfo,
                                   bool tumor_normal_mode, AltPresence& alt_presence) {
  if (support->TotalAltCov() > 0) {
    mHasAltSupport = true;
    if (tumor_normal_mode) {
      if (sinfo.TagKind() == cbdg::Label::NORMAL) alt_presence.in_normal = true;
      if (sinfo.TagKind() == cbdg::Label::TUMOR) alt_presence.in_tumor = true;
    }
  }
}

auto VariantCall::BuildPerAlleleMetrics(const VariantSupport* support, usize num_alleles) const -> PerAlleleMetrics {
  PerAlleleMetrics metrics;
  for (usize allele = 0; allele < num_alleles; ++allele) {
    const auto idx = static_cast<AlleleIndex>(allele);
    if (allele > 0) {
      metrics.ad += ","; metrics.adf += ","; metrics.adr += ","; metrics.rmq += ","; metrics.npbq += ",";
    }
    absl::StrAppend(&metrics.ad, support->TotalAlleleCov(idx));
    absl::StrAppend(&metrics.adf, support->FwdCount(idx));
    absl::StrAppend(&metrics.adr, support->RevCount(idx));
    // NOTE: absl::StrAppendFormat uses printf-style format specifiers, not fmt-style
    absl::StrAppendFormat(&metrics.rmq, "%.1f", support->RmsMappingQual(idx));
    
    // NPBQ: raw posterior base quality divided by allele depth
    // Recovers the effective per-read quality (~30 for Q30 reads at any depth)
    const auto raw_pbq = support->RawPosteriorBaseQual(idx);
    const auto allele_cov = support->TotalAlleleCov(idx);
    const auto npbq = allele_cov > 0 ? raw_pbq / static_cast<f64>(allele_cov) : 0.0;
    // NOTE: absl::StrAppendFormat uses printf-style format specifiers, not fmt-style
    absl::StrAppendFormat(&metrics.npbq, "%.1f", npbq);
  }
  return metrics;
}

// ============================================================================
// ComputeState: classify variant as SHARED, NORMAL, TUMOR, UNKNOWN, or NONE.
//
// In tumor-normal mode: SHARED/NORMAL/TUMOR based on ALT presence, NONE if no support.
// In normal-only mode: always UNKNOWN (not enough info to classify).
// ============================================================================
void VariantCall::ComputeState(const AltPresence alt_presence, const bool tumor_normal_mode) {
  static constexpr RawVariant::State STATE_MAP[4] = {
      RawVariant::State::NONE,   // 00: Neither
      RawVariant::State::NORMAL, // 01: Normal only
      RawVariant::State::TUMOR,  // 10: Tumor only
      RawVariant::State::SHARED  // 11: Both
  };

  // We compute a 2-bit mapping index directly from the boolean parameters:
  // Bit 1 (Left bit):  in_tumor
  // Bit 0 (Right bit): in_normal 
  // Evaluates exactly into [0=NONE, 1=NORMAL, 2=TUMOR, 3=SHARED].
  const usize state_idx = (static_cast<usize>(alt_presence.in_tumor) << 1) | 
                           static_cast<usize>(alt_presence.in_normal);

  mState = tumor_normal_mode ? STATE_MAP[state_idx] : RawVariant::State::UNKNOWN;
}

// ============================================================================
// BuildInfoField: assemble the VCF INFO string.
//
// Structure:  [STATE;]TYPE=<type>;LENGTH=<len>[;GRAPH_CX=...][;SEQ_CX=...]
//
//   STATE     — SHARED/NORMAL/TUMOR (tumor-normal mode only; omitted otherwise)
//   TYPE      — SNV, INS, DEL, MNP (always present)
//   LENGTH    — variant length in bp (always present)
//   GRAPH_CX  — optional graph complexity (GEI, TipToPathCovRatio, MaxDegree)
//   SEQ_CX    — optional sequence complexity (11 ML-ready features)
//
// Note: SCA, FLD, and MQCD are per-sample FORMAT fields, not site-level INFO.
// ============================================================================
void VariantCall::BuildInfoField(const bool tumor_normal_mode, FeatureFlags features) {
  using namespace std::string_view_literals;

  static constexpr std::string_view TYPE_MAP[] = {"REF"sv, "SNV"sv, "INS"sv, "DEL"sv, "MNP"sv, "CPX"sv};

  std::vector<std::string_view> vcategories;
  vcategories.reserve(mCategories.size());
  for (const auto cat : mCategories) {
    vcategories.push_back(TYPE_MAP[static_cast<i8>(cat) + 1]);
  }

  std::string info;
  info.reserve(1024);

  if (tumor_normal_mode) {
    // State prefix — tumor-normal mode only
    // SHARED/NORMAL/TUMOR state — only in tumor-normal (somatic) mode
    static constexpr std::string_view STATE_MAP[] = {"NONE"sv, "SHARED"sv, "NORMAL"sv, "TUMOR"sv, "UNKNOWN"sv};
    absl::StrAppend(&info, STATE_MAP[static_cast<i8>(mState) + 1], ";");
  }

  if (mIsMultiallelic) absl::StrAppend(&info, "MULTIALLELIC;");
  
  absl::StrAppend(&info, "TYPE=", absl::StrJoin(vcategories, ","), ";LENGTH=", absl::StrJoin(mVariantLengths, ","));
  // Optional complexity annotations — each struct owns its own VCF formatting
  if (features.enable_graph_complexity) absl::StrAppend(&info, ";GRAPH_CX=", mGraphCx.FormatVcfValue());
  if (features.enable_sequence_complexity) absl::StrAppend(&info, ";SEQ_CX=", mSeqCx.FormatVcfValue());

  mInfoField = std::move(info);
}

// ============================================================================
// AsVcfRecord: emit a VCF record with comma-separated ALTs for multi-allelic.
// ============================================================================
auto VariantCall::AsVcfRecord() const -> std::string {
  const auto alt_field = absl::StrJoin(mAltAlleles, ",");
  return fmt::format("{CHROM}\t{POS}\t.\t{REF}\t{ALT}\t{QUAL:.2f}\t.\t{INFO}\t{FORMAT}",
                     fmt::arg("CHROM", mChromName),
                     fmt::arg("POS", mStartPos1),
                     fmt::arg("REF", mRefAllele),
                     fmt::arg("ALT", alt_field),
                     fmt::arg("QUAL", mSiteQuality),
                     fmt::arg("INFO", mInfoField),
                     fmt::arg("FORMAT", absl::StrJoin(mFormatFields, "\t")));
}

}  // namespace lancet::caller
