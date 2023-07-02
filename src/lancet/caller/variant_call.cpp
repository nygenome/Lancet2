#include "lancet/caller/variant_call.h"

#include <algorithm>
#include <cmath>

#include "absl/container/btree_set.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "lancet/base/compute_stats.h"
#include "lancet/base/hash.h"
#include "lancet/hts/fisher_exact.h"
#include "spdlog/fmt/fmt.h"

namespace {

[[nodiscard]] inline auto HashRawVariant(const lancet::caller::RawVariant *var) -> u64 {
  return HashStr64(fmt::format("{},{},{},{},{},{}", var->mChromName, var->mGenomeStart1, var->mRefAllele,
                               var->mAltAllele, var->mAlleleLength, static_cast<i8>(var->mType)));
}

}  // namespace

namespace lancet::caller {

VariantCall::VariantCall(const RawVariant *var, Supports &&supprts, Samples samps, const Params &prms, const usize klen)
    : mVariantId(HashRawVariant(var)), mChromIndex(var->mChromIndex), mStartPos1(var->mGenomeStart1),
      mTotalSampleCov(0), mChromName(var->mChromName), mRefAllele(var->mRefAllele), mAltAllele(var->mAltAllele),
      mVarLength(var->mAlleleLength), mSiteQuality(0), mCategory(var->mType) {
  PerSampleEvidence per_sample_evidence;
  per_sample_evidence.reserve(supprts.size());

  for (const auto &sinfo : samps) {
    auto itr = supprts.find(sinfo.SampleName());
    if (itr == supprts.end()) {
      per_sample_evidence.emplace(sinfo, std::make_unique<VariantSupport>());
      continue;
    }

    auto handle = supprts.extract(itr);
    per_sample_evidence.emplace(sinfo, std::move(handle.mapped()));
  }

  mFormatFields.reserve(samps.size() + 1);
  mFormatFields.emplace_back("GT:AD:ADF:ADR:DP:AAF:SBS:SFS:SOS:FT:HQ:GQ:PL");

  static const auto is_normal = [](const auto &sinfo) -> bool { return sinfo.TagKind() == cbdg::Label::NORMAL; };
  static const auto is_tumor = [](const auto &sinfo) -> bool { return sinfo.TagKind() == cbdg::Label::TUMOR; };
  const auto is_germline_mode = std::ranges::all_of(samps, is_normal);
  const auto is_somatic_mode = std::ranges::any_of(samps, is_tumor);

  bool seen_in_normal = false;
  bool seen_in_tumor = false;
  std::vector<std::string> per_sample_filters;
  absl::btree_set<std::string> variant_level_filters;
  const auto is_str = var->mStrResult.mFoundStr;

  for (const auto &sinfo : samps) {
    const auto &evidence = per_sample_evidence.at(sinfo.SampleName());

    const auto sample_likelihoods = evidence->ComputePLs();
    const auto [smallest_index, second_smallest_index] = FirstAndSecondSmallestIndices(sample_likelihoods);
    const auto genotype = POSSIBLE_GENOTYPES.at(smallest_index);
    const auto gt_quality = sample_likelihoods.at(second_smallest_index);
    const auto [ref_hom_pl, het_alt_pl, alt_hom_pl] = sample_likelihoods;
    const auto [mean_ref_qual, mean_alt_qual] = evidence->MeanHaplotypeQualities();

    const auto total_alt_cov = evidence->TotalAltCov();
    const auto single_strand_alt = evidence->AltFwdCount() <= 1 || evidence->AltRevCount() <= 1;
    const auto alt_freq = evidence->AltFrequency();
    const auto strand_bias = evidence->StrandBiasScore();
    const auto fet_score = SomaticFetScore(sinfo, per_sample_evidence);
    const auto odds_score = SomaticOddsScore(sinfo, per_sample_evidence);

    const auto somatic_quality = static_cast<int>(std::min(fet_score, odds_score));
    mSiteQuality = std::max(mSiteQuality, is_germline_mode ? gt_quality : somatic_quality);
    mTotalSampleCov += evidence->TotalSampleCov();
    per_sample_filters.clear();

    if (sinfo.TagKind() == cbdg::Label::NORMAL) {
      // NOLINTBEGIN(readability-braces-around-statements)
      if (mTotalSampleCov < prms.mMinNmlCov) per_sample_filters.emplace_back("LowNmlCov");
      // NOLINTEND(readability-braces-around-statements)
      if (is_germline_mode && genotype != REF_HOM && (single_strand_alt || strand_bias > prms.mMinPhredScore)) {
        per_sample_filters.emplace_back("StrandBias");
      }
    }

    if (sinfo.TagKind() == cbdg::Label::TUMOR) {
      // NOLINTBEGIN(readability-braces-around-statements)
      if (mTotalSampleCov < prms.mMinTmrCov) per_sample_filters.emplace_back("LowTmrCov");
      if (alt_freq < prms.mMinTmrVaf) per_sample_filters.emplace_back("LowTmrVaf");
      if (total_alt_cov < prms.mMinTmrAltCnt) per_sample_filters.emplace_back("LowTmrCnt");
      if (fet_score < prms.mMinPhredScore) per_sample_filters.emplace_back("LowSomaticFS");
      if (odds_score < prms.mMinPhredScore) per_sample_filters.emplace_back("LowSomaticOdds");
      // NOLINTEND(readability-braces-around-statements)
      if (genotype != REF_HOM && (single_strand_alt || strand_bias > prms.mMinPhredScore)) {
        per_sample_filters.emplace_back("StrandBias");
      }
    }

    std::ranges::sort(per_sample_filters);
    variant_level_filters.insert(per_sample_filters.cbegin(), per_sample_filters.cend());
    const auto sample_ft = per_sample_filters.empty() ? "PASS" : absl::StrJoin(per_sample_filters, ";");

    // NOLINTBEGIN(readability-braces-around-statements)
    if (genotype != REF_HOM && sinfo.TagKind() == cbdg::Label::NORMAL) seen_in_normal = true;
    if (genotype != REF_HOM && sinfo.TagKind() == cbdg::Label::TUMOR) seen_in_tumor = true;
    // NOLINTEND(readability-braces-around-statements)

    mFormatFields.emplace_back(
        // GT:AD:ADF:ADR:DP:VAF:SBS:SFS:SOS:FT:HQ:GQ:PL
        fmt::format("{}:{},{}:{},{}:{},{}:{}:{:.4f}:{}:{}:{}:{}:{},{}:{}:{},{},{}", genotype, evidence->TotalRefCov(),
                    evidence->TotalAltCov(), evidence->RefFwdCount(), evidence->AltFwdCount(), evidence->RefRevCount(),
                    evidence->AltRevCount(), evidence->TotalSampleCov(), alt_freq, strand_bias, fet_score, odds_score,
                    sample_ft, mean_ref_qual, mean_alt_qual, gt_quality, ref_hom_pl, het_alt_pl, alt_hom_pl));
  }

  mFilterField = variant_level_filters.empty() ? "PASS" : absl::StrJoin(variant_level_filters, ";");
  mState = is_somatic_mode && seen_in_tumor && mFilterField == "PASS" ? RawVariant::State::SOMATIC
           : seen_in_normal && seen_in_tumor                          ? RawVariant::State::SHARED
           : seen_in_normal                                           ? RawVariant::State::NORMAL
                                                                      : RawVariant::State::NONE;

  using namespace std::string_view_literals;
  const auto vstate = mState == RawVariant::State::SHARED    ? "SHARED"sv
                      : mState == RawVariant::State::NORMAL  ? "NORMAL"sv
                      : mState == RawVariant::State::SOMATIC ? "SOMATIC"sv
                                                             : "NONE"sv;

  const auto vcategory = mCategory == RawVariant::Type::SNV   ? "SNV"sv
                         : mCategory == RawVariant::Type::INS ? "INS"sv
                         : mCategory == RawVariant::Type::DEL ? "DEL"sv
                         : mCategory == RawVariant::Type::MNP ? "MNP"sv
                                                              : "REF"sv;

  const auto str_info = is_str ? fmt::format(";STR={}{}", var->mStrResult.mStrLen, var->mStrResult.mStrMotif) : "";
  mInfoField = fmt::format("{};CATEGORY={};LEN={};KMERSIZE={}{}", vstate, vcategory, mVarLength, klen, str_info);
}

auto VariantCall::AsVcfRecord() const -> std::string {
  auto vcf_record = fmt::format("{}\t{}\t.\t{}\t{}\t{}\t{}\t{}\t{}", mChromName, mStartPos1, mRefAllele, mAltAllele,
                                mSiteQuality, mFilterField, mInfoField, absl::StrJoin(mFormatFields, "\t"));

  // No newline. caller of this method will add new line if needed
  return vcf_record;
}

auto VariantCall::SomaticFetScore(const core::SampleInfo &current, const PerSampleEvidence &supports) -> u8 {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (current.TagKind() != cbdg::Label::TUMOR) return 0;

  u32 current_tmr_alt = 0;
  u32 current_tmr_ref = 0;
  OnlineStats nml_alts;
  OnlineStats nml_refs;

  for (const auto &[sample_info, evidence] : supports) {
    // Ignore other tumor samples even if present
    if (sample_info.SampleName() == current.SampleName()) {
      current_tmr_alt = evidence->TotalAltCov();
      current_tmr_ref = evidence->TotalRefCov();
      continue;
    }

    if (sample_info.TagKind() == cbdg::Label::NORMAL) {
      nml_alts.Add(evidence->TotalAltCov());
      nml_refs.Add(evidence->TotalRefCov());
    }
  }

  const auto cnt_tmr_alt = static_cast<int>(current_tmr_alt);
  const auto cnt_tmr_ref = static_cast<int>(current_tmr_ref);
  const auto avg_nml_alt = static_cast<int>(std::round(nml_alts.Mean()));
  const auto avg_nml_ref = static_cast<int>(std::round(nml_refs.Mean()));

  using Row = hts::FisherExact::Row;
  const auto tmr_counts = Row{cnt_tmr_alt, cnt_tmr_ref};
  const auto nml_counts = Row{avg_nml_alt, avg_nml_ref};
  const auto result = hts::FisherExact::Test({tmr_counts, nml_counts});
  return hts::ErrorProbToPhred(result.mMoreProb);
}

auto VariantCall::SomaticOddsScore(const core::SampleInfo &current, const PerSampleEvidence &supports) -> u8 {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (current.TagKind() != cbdg::Label::TUMOR) return 0;

  f64 nml_vaf = std::numeric_limits<f64>::min();
  f64 tmr_vaf = std::numeric_limits<f64>::min();

  for (const auto &[sample_info, evidence] : supports) {
    // Ignore other tumor samples even if present
    if (sample_info.SampleName() == current.SampleName()) {
      tmr_vaf = evidence->AltFrequency();
      continue;
    }

    if (sample_info.TagKind() == cbdg::Label::NORMAL) {
      nml_vaf = std::max(nml_vaf, evidence->AltFrequency());
    }
  }

  const auto odds_ratio = tmr_vaf / nml_vaf;
  const auto somatic_error_prob = 1.0 / (odds_ratio + 1.0);
  return hts::ErrorProbToPhred(somatic_error_prob);
}

auto VariantCall::FirstAndSecondSmallestIndices(const std::array<int, 3> &pls) -> std::array<usize, 2> {
  const auto *min_itr = std::ranges::min_element(pls);
  const auto smallest_idx = static_cast<usize>(std::distance(pls.cbegin(), min_itr));

  const auto second_min_comparator = [&min_itr](const int lhs, const int rhs) -> bool {
    return (lhs == *min_itr) ? false : (rhs == *min_itr) ? true : lhs < rhs;
  };

  const auto *second_itr = std::ranges::min_element(pls, second_min_comparator);
  const auto second_smallest_idx = static_cast<usize>(std::distance(pls.cbegin(), second_itr));

  return {smallest_idx, second_smallest_idx};
}

}  // namespace lancet::caller
