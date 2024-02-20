#include "lancet/caller/variant_call.h"

#include <algorithm>
#include <cmath>

#include "absl/container/btree_set.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "lancet/base/compute_stats.h"
#include "lancet/base/hash.h"
#include "lancet/hts/fisher_exact.h"
#include "lancet/hts/phred_quality.h"
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
  mFormatFields.emplace_back("GT:AD:ADF:ADR:DP:WDC:WTC:PRF:VAF:AFR:SFS:AQM:AQR:MQM:MQR:ASDM:ASDR:GQ:PL");

  static const auto is_normal = [](const auto &sinfo) -> bool { return sinfo.TagKind() == cbdg::Label::NORMAL; };
  const auto germline_mode = std::ranges::all_of(samps, is_normal);

  bool alt_seen_in_normal = false;
  bool alt_seen_in_tumor = false;
  std::vector<std::string> current_filters;
  absl::btree_set<std::string> variant_site_filters;
  const auto is_str = var->mStrResult.mFoundStr;

  for (const auto &sinfo : samps) {
    const auto &evidence = per_sample_evidence.at(sinfo.SampleName());

    const auto phred_likelihoods = evidence->ComputePLs();
    const auto [ref_hom_pl, het_alt_pl, alt_hom_pl] = phred_likelihoods;
    const auto [smallest_index, second_smallest_index] = FirstAndSecondSmallestIndices(phred_likelihoods);

    const auto genotype = POSSIBLE_GENOTYPES.at(smallest_index);
    const auto genotype_quality = static_cast<u32>(phred_likelihoods.at(second_smallest_index));

    const auto allele_qual_stats = evidence->AlleleQualityStats();
    const auto [raq_median, aaq_median] = allele_qual_stats.ref_alt_medians;
    const auto [raq_range, aaq_range] = allele_qual_stats.ref_alt_ranges;

    const auto mapping_qual_stats = evidence->MappingQualityStats();
    const auto [rmq_median, amq_median] = mapping_qual_stats.ref_alt_medians;
    const auto [rmq_range, amq_range] = mapping_qual_stats.ref_alt_ranges;

    const auto aln_score_stats = evidence->AlnDiffScoreStats();
    const auto [rasd_median, aasd_median] = aln_score_stats.ref_alt_medians;
    const auto [rasd_range, aasd_range] = aln_score_stats.ref_alt_ranges;

    const auto alt_allele_freq = evidence->AltFrequency();
    const auto alt_on_single_strand = evidence->AltFwdCount() == 0 || evidence->AltRevCount() == 0;
    const auto tn_alt_odds_ratio = SomaticOddsRatio(sinfo, per_sample_evidence);
    const auto fisher_score = SomaticFisherScore(sinfo, per_sample_evidence);

    mSiteQuality = std::max(mSiteQuality, germline_mode ? static_cast<f64>(ref_hom_pl) : fisher_score);
    mTotalSampleCov += evidence->TotalSampleCov();
    current_filters.clear();

    if (sinfo.TagKind() == cbdg::Label::NORMAL) {
      // NOLINTBEGIN(readability-braces-around-statements)
      if (evidence->TotalSampleCov() < prms.mMinNmlCov) current_filters.emplace_back("LowNmlCov");
      if (germline_mode && genotype != REF_HOM && alt_on_single_strand) current_filters.emplace_back("StrandBias");
      // NOLINTEND(readability-braces-around-statements)
    }

    if (sinfo.TagKind() == cbdg::Label::TUMOR) {
      using RawVariant::Type::SNV;
      // NOLINTBEGIN(readability-braces-around-statements)
      if (evidence->TotalSampleCov() < prms.mMinTmrCov) current_filters.emplace_back("LowTmrCov");
      if (genotype != REF_HOM && alt_on_single_strand) current_filters.emplace_back("StrandBias");
      if (tn_alt_odds_ratio < prms.mMinOddsRatio) current_filters.emplace_back("LowOddsRatio");
      if (is_str && fisher_score < prms.mMinStrFisher) current_filters.emplace_back("LowStrFisher");
      if (mCategory != SNV && fisher_score < prms.mMinInDelFisher) current_filters.emplace_back("LowIndelFisher");
      if (mCategory == SNV && fisher_score < prms.mMinSnvFisher) current_filters.emplace_back("LowSnvFisher");
      // NOLINTEND(readability-braces-around-statements)
    }

    std::ranges::sort(current_filters);
    variant_site_filters.insert(current_filters.cbegin(), current_filters.cend());

    // NOLINTBEGIN(readability-braces-around-statements)
    if (genotype != REF_HOM && sinfo.TagKind() == cbdg::Label::NORMAL) alt_seen_in_normal = true;
    if (genotype != REF_HOM && sinfo.TagKind() == cbdg::Label::TUMOR) alt_seen_in_tumor = true;
    // NOLINTEND(readability-braces-around-statements)

    mFormatFields.emplace_back(fmt::format(
        "{GT}:{REF_AD},{ALT_AD}:{REF_ADF},{ALT_ADF}:{REF_ADR},{ALT_ADR}:{DP}:{WDC:.2f}:{WTC:.2f}:{PRF:.2f}:"
        "{VAF:.2f}:{AFR:.2f}:{SFS:.2f}:{REF_AQM},{ALT_AQM}:{REF_AQR},{ALT_AQR}:{REF_MQM},{ALT_MQM}:{REF_"
        "MQR},{ALT_MQR}:{REF_ASDM},{ALT_ASDM}:{REF_ASDR},{ALT_ASDR}:{GQ}:{PL1},{PL2},{PL3}",
        fmt::arg("GT", genotype), fmt::arg("REF_AD", evidence->TotalRefCov()),
        fmt::arg("ALT_AD", evidence->TotalAltCov()), fmt::arg("REF_ADF", evidence->RefFwdCount()),
        fmt::arg("ALT_ADF", evidence->AltFwdCount()), fmt::arg("REF_ADR", evidence->RefRevCount()),
        fmt::arg("ALT_ADR", evidence->AltRevCount()), fmt::arg("DP", evidence->TotalSampleCov()),
        fmt::arg("WDC", sinfo.MeanSampledCov()), fmt::arg("WTC", sinfo.MeanTotalCov()),
        fmt::arg("PRF", sinfo.PassReadsFraction()), fmt::arg("VAF", alt_allele_freq),
        fmt::arg("AFR", tn_alt_odds_ratio), fmt::arg("SFS", fisher_score), fmt::arg("REF_AQM", raq_median),
        fmt::arg("ALT_AQM", aaq_median), fmt::arg("REF_AQR", raq_range), fmt::arg("ALT_AQR", aaq_range),
        fmt::arg("REF_MQM", rmq_median), fmt::arg("ALT_MQM", amq_median), fmt::arg("REF_MQR", rmq_range),
        fmt::arg("ALT_MQR", amq_range), fmt::arg("REF_ASDM", rasd_median), fmt::arg("ALT_ASDM", aasd_median),
        fmt::arg("REF_ASDR", rasd_range), fmt::arg("ALT_ASDR", aasd_range), fmt::arg("GQ", genotype_quality),
        fmt::arg("PL1", ref_hom_pl), fmt::arg("PL2", het_alt_pl), fmt::arg("PL3", alt_hom_pl)));
  }

  mFilterField = variant_site_filters.empty() ? "PASS" : absl::StrJoin(variant_site_filters, ";");
  mState = alt_seen_in_normal && alt_seen_in_tumor ? RawVariant::State::SHARED
           : alt_seen_in_normal                    ? RawVariant::State::NORMAL
           : alt_seen_in_tumor                     ? RawVariant::State::TUMOR
                                                   : RawVariant::State::NONE;

  using namespace std::string_view_literals;
  const auto vstate = mState == RawVariant::State::SHARED   ? "SHARED"sv
                      : mState == RawVariant::State::NORMAL ? "NORMAL"sv
                      : mState == RawVariant::State::TUMOR  ? "TUMOR"sv
                                                            : "NONE"sv;

  const auto vcategory = mCategory == RawVariant::Type::SNV   ? "SNV"sv
                         : mCategory == RawVariant::Type::INS ? "INS"sv
                         : mCategory == RawVariant::Type::DEL ? "DEL"sv
                         : mCategory == RawVariant::Type::MNP ? "MNP"sv
                                                              : "REF"sv;

  const auto str_info = is_str ? fmt::format(";STR={}{}", var->mStrResult.mStrLen, var->mStrResult.mStrMotif) : "";
  mInfoField = fmt::format("{};TYPE={};LEN={};KMER_LEN={}{}", vstate, vcategory, mVarLength, klen, str_info);
}

auto VariantCall::AsVcfRecord() const -> std::string {
  auto vcf_record = fmt::format("{}\t{}\t.\t{}\t{}\t{:.4f}\t{}\t{}\t{}", mChromName, mStartPos1, mRefAllele, mAltAllele,
                                mSiteQuality, mFilterField, mInfoField, absl::StrJoin(mFormatFields, "\t"));

  // No newline. caller of this method will add new line if needed
  return vcf_record;
}

auto VariantCall::SomaticFisherScore(const core::SampleInfo &curr, const PerSampleEvidence &supports) -> f64 {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (curr.TagKind() != cbdg::Label::TUMOR) return 0;

  u32 current_tmr_alt = 0;
  u32 current_tmr_ref = 0;
  OnlineStats nml_alts;
  OnlineStats nml_refs;

  for (const auto &[sample_info, evidence] : supports) {
    // Ignore other tumor samples even if present
    if (sample_info.SampleName() == curr.SampleName()) {
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

auto VariantCall::SomaticOddsRatio(const core::SampleInfo &curr, const PerSampleEvidence &supports) -> f64 {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (curr.TagKind() != cbdg::Label::TUMOR) return 0;

  f64 max_nml_vaf = std::numeric_limits<f64>::min();
  f64 curr_tmr_vaf = std::numeric_limits<f64>::min();

  for (const auto &[sample_info, evidence] : supports) {
    // Ignore other tumor samples even if present
    if (sample_info.SampleName() == curr.SampleName()) {
      curr_tmr_vaf = evidence->AltFrequency();
      continue;
    }

    if (sample_info.TagKind() == cbdg::Label::NORMAL) {
      max_nml_vaf = std::max(max_nml_vaf, evidence->AltFrequency());
    }
  }

  const auto odds_ratio = curr_tmr_vaf / max_nml_vaf;
  return std::clamp(odds_ratio, 0.0, static_cast<f64>(hts::MAX_PHRED_SCORE));
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
