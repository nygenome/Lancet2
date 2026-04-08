#include "lancet/caller/variant_call.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "lancet/base/assert.h"
#include "lancet/base/compute_stats.h"
#include "lancet/base/hash.h"
#include "lancet/base/longdust_scorer.h"
#include "lancet/base/polar_coords.h"
#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_support.h"
#include "lancet/cbdg/label.h"
#include "spdlog/fmt/bundled/core.h"

namespace {

[[nodiscard]] inline auto HashRawVariant(const lancet::caller::RawVariant* var) -> u64 {
  return HashStr64(fmt::format("{},{},{},{},{},{}", var->mChromName, var->mGenomeStart1, var->mRefAllele,
                               var->mAltAllele, var->mAlleleLength, static_cast<i8>(var->mType)));
}

[[nodiscard]] inline auto HashVariantGroup(absl::Span<const lancet::caller::RawVariant* const> vars) -> u64 {
  std::string combined;
  combined.reserve(vars.size() * 80);
  for (const auto* var : vars) {
    absl::StrAppend(&combined, var->mChromName, ",", var->mGenomeStart1, ",", var->mRefAllele, ",",
                     var->mAltAllele, ",", var->mAlleleLength, ",", static_cast<i8>(var->mType), ",");
  }
  return HashStr64(combined);
}

}  // namespace

namespace lancet::caller {

// ============================================================================
// Multi-allelic constructor
//
// Takes a group of RawVariants at the same locus and merges their per-sample
// evidence. Each RawVariant's AlleleIndex mapping (from genotyper) is used
// to build a unified multi-allelic VariantSupport per sample.
// ============================================================================
VariantCall::VariantCall(VariantGroup variants, SupportsByVariant&& all_supports, Samples samps,
                         FeatureFlags features, const f64 window_cov)
    : mVariantId(HashVariantGroup(variants)),
      mChromIndex(variants[0]->mChromIndex), mStartPos1(variants[0]->mGenomeStart1),
      mTotalSampleCov(0), mChromName(variants[0]->mChromName), mRefAllele(variants[0]->mRefAllele),
      mVariantLength(variants[0]->mAlleleLength), mSiteQuality(0), mCategory(variants[0]->mType),
      mFeatureFlags(features), mWindowCov(window_cov) {
  // Collect ALT alleles and merge complexity scores (element-wise max across the group)
  for (const auto* var : variants) {
    mAltAlleles.push_back(var->mAltAllele);

    // Graph metrics: take from first variant (all share same window-level metrics)
    if (var == variants[0]) {
      mGraphMetrics = var->mGraphMetrics;
    }

    // Sequence complexity: element-wise max across all variants at this locus
    mSeqCx.MergeMax(var->mSeqCx);
  }

  // Merge per-variant bi-allelic evidence into a single multi-allelic VariantSupport.
  //
  // Each variant independently has alleles {0=REF, 1=ALT}. When merging N
  // variants at a locus, we remap each one's ALT to its position in the
  // multi-allelic allele list:
  //
  //   variants[0] ALT(1) → merged allele 1
  //   variants[1] ALT(1) → merged allele 2
  //   variants[2] ALT(1) → merged allele 3 ...
  //
  // REF reads are shared — we merge REF(0) from all variants into merged allele 0
  // (deduplication by read name hash prevents double-counting).
  SupportArray per_sample_evidence;

  for (const auto& sinfo : samps) {
    auto merged = std::make_unique<VariantSupport>();

    for (usize var_idx = 0; var_idx < variants.size(); ++var_idx) {
      auto var_it = all_supports.find(variants[var_idx]);
      if (var_it == all_supports.end()) continue;

      const auto* src = var_it->second.Find(sinfo.SampleName());
      if (src == nullptr) continue;

      // Merge REF reads (allele 0 → 0) from all variants, dedup handles overlaps
      merged->MergeAlleleFrom(*src, REF_ALLELE_IDX, REF_ALLELE_IDX);
      // Remap this variant's ALT (allele 1) → its multi-allelic position (var_idx + 1)
      merged->MergeAlleleFrom(*src, AlleleIndex{1}, static_cast<AlleleIndex>(var_idx + 1));
    }

    per_sample_evidence.Insert(sinfo.SampleName(), std::move(merged));
  }

  Finalize(per_sample_evidence, samps, features);
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

    const auto pls = support->ComputePLs();
    const auto gq_value = VariantSupport::ComputeGQ(pls);

    // Best genotype = minimum PL index
    usize best_gt_idx = 0;
    if (!pls.empty()) {
      best_gt_idx = static_cast<usize>(
          std::distance(pls.cbegin(), std::min_element(pls.cbegin(), pls.cend())));
    }
    const auto genotype = GenotypeFromGLIndex(best_gt_idx, num_alleles);

    // Site quality: somatic log odds ratio in tumor-normal mode,
    // per-read evidence strength (PL[0]/DP, capped at 10) in normal-only mode.
    const auto somatic_lor = SomaticLogOddsRatio(sinfo, evidence, samps);
    const auto ref_hom_pl = pls.empty() ? 0 : pls[0];
    const auto sample_dp = support->TotalSampleCov();
    // Per-read evidence against hom-REF: PL[0]/DP, capped at 10.0.
    // Coverage-invariant: clean het returns ~3.0 at any depth.
    // Take max across multiple normals (conservative: strongest evidence wins).
    const auto per_read_qual = (!pls.empty() && sample_dp > 0)
        ? std::min(static_cast<f64>(ref_hom_pl) / static_cast<f64>(sample_dp), 10.0)
        : 0.0;
    mSiteQuality = std::max(mSiteQuality, tumor_normal_mode ? somatic_lor : per_read_qual);
    mTotalSampleCov += support->TotalSampleCov();

    // Track ALT support globally (for HasAltSupport()) and per-label (for state)
    if (support->TotalAltCov() > 0) {
      mHasAltSupport = true;
      if (tumor_normal_mode) {
        // NOLINTBEGIN(readability-braces-around-statements)
        if (sinfo.TagKind() == cbdg::Label::NORMAL) alt_presence.in_normal = true;
        if (sinfo.TagKind() == cbdg::Label::TUMOR) alt_presence.in_tumor = true;
        // NOLINTEND(readability-braces-around-statements)
      }
    }

    // Per-allele metric strings (Number=R fields)
    std::vector<std::string> ad_vals, adf_vals, adr_vals, rmq_vals, npbq_vals;
    ad_vals.reserve(num_alleles);
    adf_vals.reserve(num_alleles);
    adr_vals.reserve(num_alleles);
    rmq_vals.reserve(num_alleles);
    npbq_vals.reserve(num_alleles);

    for (usize allele = 0; allele < num_alleles; ++allele) {
      const auto idx = static_cast<AlleleIndex>(allele);
      ad_vals.push_back(std::to_string(support->TotalAlleleCov(idx)));
      adf_vals.push_back(std::to_string(support->FwdCount(idx)));
      adr_vals.push_back(std::to_string(support->RevCount(idx)));
      rmq_vals.push_back(fmt::format("{:.1f}", support->RmsMappingQual(idx)));
      // NPBQ: raw posterior base quality divided by allele depth
      // Recovers the effective per-read quality (~30 for Q30 reads at any depth)
      const auto raw_pbq = support->RawPosteriorBaseQual(idx);
      const auto allele_cov = support->TotalAlleleCov(idx);
      const auto npbq = allele_cov > 0 ? raw_pbq / static_cast<f64>(allele_cov) : 0.0;
      npbq_vals.push_back(fmt::format("{:.1f}", npbq));
    }

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

    // PL string (Number=G field)
    std::vector<std::string> pl_strs;
    pl_strs.reserve(pls.size());
    for (const auto pl_val : pls) {
      pl_strs.push_back(std::to_string(pl_val));
    }

    mFormatFields.emplace_back(fmt::format(
        "{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
        genotype,
        absl::StrJoin(ad_vals, ","),
        absl::StrJoin(adf_vals, ","),
        absl::StrJoin(adr_vals, ","),
        support->TotalSampleCov(),
        absl::StrJoin(rmq_vals, ","),
        absl::StrJoin(npbq_vals, ","),
        sb_val,
        sca_val,
        fld_val,
        rpcd_val,
        bqcd_val,
        mqcd_val,
        asmd_val,
        sdfc_val,
        prad_val,
        pang_val,
        pls.empty() ? "." : absl::StrJoin(pl_strs, ","),
        gq_value));
  }

  return alt_presence;
}

// ============================================================================
// ComputeState: classify variant as SHARED, NORMAL, TUMOR, UNKNOWN, or NONE.
//
// In tumor-normal mode: SHARED/NORMAL/TUMOR based on ALT presence, NONE if no support.
// In normal-only mode: always UNKNOWN (not enough info to classify).
// ============================================================================
void VariantCall::ComputeState(const AltPresence alt_presence, const bool tumor_normal_mode) {
  if (!tumor_normal_mode) {
    mState = RawVariant::State::UNKNOWN;
    return;
  }

  // NOLINTBEGIN(readability-avoid-nested-conditional-operator)
  mState = alt_presence.in_normal && alt_presence.in_tumor ? RawVariant::State::SHARED
           : alt_presence.in_normal                        ? RawVariant::State::NORMAL
           : alt_presence.in_tumor                         ? RawVariant::State::TUMOR
                                                            : RawVariant::State::NONE;
  // NOLINTEND(readability-avoid-nested-conditional-operator)
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

  // NOLINTBEGIN(readability-avoid-nested-conditional-operator)
  const auto vcategory = mCategory == RawVariant::Type::SNV   ? "SNV"sv
                         : mCategory == RawVariant::Type::INS ? "INS"sv
                         : mCategory == RawVariant::Type::DEL ? "DEL"sv
                         : mCategory == RawVariant::Type::MNP ? "MNP"sv
                                                              : "REF"sv;
  // NOLINTEND(readability-avoid-nested-conditional-operator)

  std::string info;
  info.reserve(1024);

  // State prefix — tumor-normal mode only
  if (tumor_normal_mode) {
    // SHARED/NORMAL/TUMOR state — only in tumor-normal (somatic) mode
    // NOLINTBEGIN(readability-avoid-nested-conditional-operator)
    const auto vstate = mState == RawVariant::State::SHARED ? "SHARED"sv
                        : mState == RawVariant::State::NORMAL ? "NORMAL"sv
                        : mState == RawVariant::State::TUMOR  ? "TUMOR"sv
                                                              : "NONE"sv;
    // NOLINTEND(readability-avoid-nested-conditional-operator)
    absl::StrAppend(&info, vstate, ";");
  }

  absl::StrAppend(&info, "TYPE=", vcategory, ";LENGTH=", mVariantLength);

  // Optional complexity annotations — each struct owns its own VCF formatting
  if (features.enable_graph_complexity) {
    absl::StrAppend(&info, ";GRAPH_CX=", mGraphMetrics.FormatVcfValue());
  }

  if (features.enable_sequence_complexity) {
    absl::StrAppend(&info, ";SEQ_CX=", mSeqCx.FormatVcfValue());
  }

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
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (curr.TagKind() != cbdg::Label::TUMOR) return 0;

  u32 current_tmr_alt = 0;
  u32 current_tmr_ref = 0;
  OnlineStats nml_alts;
  OnlineStats nml_refs;

  for (const auto& sinfo : samps) {
    if (sinfo.SampleName() == curr.SampleName()) {
      const auto* evidence = supports.Find(sinfo.SampleName());
      if (evidence != nullptr) {
        current_tmr_alt = evidence->TotalAltCov();
        current_tmr_ref = evidence->TotalRefCov();
      }
      continue;
    }

    if (sinfo.TagKind() == cbdg::Label::NORMAL) {
      const auto* evidence = supports.Find(sinfo.SampleName());
      if (evidence != nullptr) {
        nml_alts.Add(evidence->TotalAltCov());
        nml_refs.Add(evidence->TotalRefCov());
      }
    }
  }

  // Haldane correction (+1 to all cells)
  const auto ta = static_cast<f64>(current_tmr_alt) + 1.0;
  const auto tr = static_cast<f64>(current_tmr_ref) + 1.0;
  const auto na = std::round(nml_alts.Mean()) + 1.0;
  const auto nr = std::round(nml_refs.Mean()) + 1.0;

  return std::log((ta * nr) / (tr * na));
}

}  // namespace lancet::caller
