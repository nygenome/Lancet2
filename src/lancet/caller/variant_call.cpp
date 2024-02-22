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

VariantCall::VariantCall(const RawVariant *var, Supports &&supports, Samples samps, const usize kmerlen)
    : mVariantId(HashRawVariant(var)), mChromIndex(var->mChromIndex), mStartPos1(var->mGenomeStart1),
      mTotalSampleCov(0), mChromName(var->mChromName), mRefAllele(var->mRefAllele), mAltAllele(var->mAltAllele),
      mVariantLength(var->mAlleleLength), mSiteQuality(0), mCategory(var->mType) {
  PerSampleEvidence per_sample_evidence;
  per_sample_evidence.reserve(supports.size());

  for (const auto &sinfo : samps) {
    auto itr = supports.find(sinfo.SampleName());
    if (itr == supports.end()) {
      per_sample_evidence.emplace(sinfo, std::make_unique<VariantSupport>());
      continue;
    }

    auto handle = supports.extract(itr);
    per_sample_evidence.emplace(sinfo, std::move(handle.mapped()));
  }

  mFormatFields.reserve(samps.size() + 1);
  mFormatFields.emplace_back("GT:AD:ADF:ADR:DP:WDC:WTC:PRF:VAF:RAQS:AAQS:RMQS:AMQS:RAPDS:AAPDS:GQ:PL");

  static const auto is_normal = [](const auto &sinfo) -> bool { return sinfo.TagKind() == cbdg::Label::NORMAL; };
  const auto germline_mode = std::ranges::all_of(samps, is_normal);

  bool alt_seen_in_normal = false;
  bool alt_seen_in_tumor = false;
  const auto is_str = var->mStrResult.mFoundStr;

  for (const auto &sinfo : samps) {
    const auto &evidence = per_sample_evidence.at(sinfo.SampleName());

    const auto phred_likelihoods = evidence->ComputePLs();
    const auto [ref_hom_pl, het_alt_pl, alt_hom_pl] = phred_likelihoods;
    const auto [smallest_index, second_smallest_index] = FirstAndSecondSmallestIndices(phred_likelihoods);

    const auto genotype = POSSIBLE_GENOTYPES.at(smallest_index);
    const auto genotype_quality = static_cast<u32>(phred_likelihoods.at(second_smallest_index));

    const auto allele_qual_stats = evidence->AlleleQualityStats();
    const auto mapping_qual_stats = evidence->MappingQualityStats();
    const auto aln_score_stats = evidence->AlnDiffScoreStats();

    const auto alt_allele_freq = evidence->AltFrequency();
    const auto fisher_score = SomaticFisherScore(sinfo, per_sample_evidence);

    mSiteQuality = std::max(mSiteQuality, germline_mode ? static_cast<f64>(ref_hom_pl) : fisher_score);
    mTotalSampleCov += evidence->TotalSampleCov();

    // NOLINTBEGIN(readability-braces-around-statements)
    if (genotype != REF_HOM && sinfo.TagKind() == cbdg::Label::NORMAL) alt_seen_in_normal = true;
    if (genotype != REF_HOM && sinfo.TagKind() == cbdg::Label::TUMOR) alt_seen_in_tumor = true;
    // NOLINTEND(readability-braces-around-statements)

    mFormatFields.emplace_back(fmt::format(
        "{GT}:{AD1},{AD2}:{ADF1},{ADF2}:{ADR1},{ADR2}:"
        "{DP}:{WDC:.2f}:{WTC:.2f}:{PRF:.2f}:{VAF:.2f}:"
        "{RAQ_MIN},{RAQ_MEDIAN},{RAQ_MAX},{RAQ_MAD}:"
        "{AAQ_MIN},{AAQ_MEDIAN},{AAQ_MAX},{AAQ_MAD}:"
        "{RMQ_MIN},{RMQ_MEDIAN},{RMQ_MAX},{RMQ_MAD}:"
        "{AMQ_MIN},{AMQ_MEDIAN},{AMQ_MAX},{AMQ_MAD}:"
        "{RAPD_MIN},{RAPD_MEDIAN},{RAPD_MAX},{RAPD_MAD}:"
        "{AAPD_MIN},{AAPD_MEDIAN},{AAPD_MAX},{AAPD_MAD}:"
        "{GQ}:{HOM_REF_PL},{HET_ALT_PL},{HOM_ALT_PL}",

        fmt::arg("GT", genotype),

        fmt::arg("AD1", evidence->TotalRefCov()), fmt::arg("AD2", evidence->TotalAltCov()),
        fmt::arg("ADF1", evidence->RefFwdCount()), fmt::arg("ADF2", evidence->AltFwdCount()),
        fmt::arg("ADR1", evidence->RefRevCount()), fmt::arg("ADR2", evidence->AltRevCount()),

        fmt::arg("DP", evidence->TotalSampleCov()), fmt::arg("WDC", sinfo.MeanSampledCov()),
        fmt::arg("WTC", sinfo.MeanTotalCov()), fmt::arg("PRF", sinfo.PassReadsFraction()),
        fmt::arg("VAF", alt_allele_freq),

        fmt::arg("RAQ_MIN", allele_qual_stats.refMinVal), fmt::arg("RAQ_MEDIAN", allele_qual_stats.refMedian),
        fmt::arg("RAQ_MAX", allele_qual_stats.refMaxVal), fmt::arg("RAQ_MAD", allele_qual_stats.refMADVal),

        fmt::arg("AAQ_MIN", allele_qual_stats.altMinVal), fmt::arg("AAQ_MEDIAN", allele_qual_stats.altMedian),
        fmt::arg("AAQ_MAX", allele_qual_stats.altMaxVal), fmt::arg("AAQ_MAD", allele_qual_stats.altMADVal),

        fmt::arg("RMQ_MIN", mapping_qual_stats.refMinVal), fmt::arg("RMQ_MEDIAN", mapping_qual_stats.refMedian),
        fmt::arg("RMQ_MAX", mapping_qual_stats.refMaxVal), fmt::arg("RMQ_MAD", mapping_qual_stats.refMADVal),

        fmt::arg("AMQ_MIN", mapping_qual_stats.altMinVal), fmt::arg("AMQ_MEDIAN", mapping_qual_stats.altMedian),
        fmt::arg("AMQ_MAX", mapping_qual_stats.altMaxVal), fmt::arg("AMQ_MAD", mapping_qual_stats.altMADVal),

        fmt::arg("RAPD_MIN", aln_score_stats.refMinVal), fmt::arg("RAPD_MEDIAN", aln_score_stats.refMedian),
        fmt::arg("RAPD_MAX", aln_score_stats.refMaxVal), fmt::arg("RAPD_MAD", aln_score_stats.refMADVal),

        fmt::arg("AAPD_MIN", aln_score_stats.altMinVal), fmt::arg("AAPD_MEDIAN", aln_score_stats.altMedian),
        fmt::arg("AAPD_MAX", aln_score_stats.altMaxVal), fmt::arg("AAPD_MAD", aln_score_stats.altMADVal),

        fmt::arg("GQ", genotype_quality),

        fmt::arg("HOM_REF_PL", ref_hom_pl), fmt::arg("HET_ALT_PL", het_alt_pl), fmt::arg("HOM_ALT_PL", alt_hom_pl)));
  }

  mState = alt_seen_in_normal && alt_seen_in_tumor ? RawVariant::State::SHARED
           : alt_seen_in_normal                    ? RawVariant::State::NORMAL
           : alt_seen_in_tumor                     ? RawVariant::State::TUMOR
                                                   : RawVariant::State::NONE;

  using namespace std::string_view_literals;
  const auto vstate = alt_seen_in_normal && alt_seen_in_tumor ? "SHARED"sv
                      : alt_seen_in_normal                    ? "NORMAL"sv
                      : alt_seen_in_tumor                     ? "TUMOR"sv
                                                              : "NONE"sv;

  const auto vcategory = mCategory == RawVariant::Type::SNV   ? "SNV"sv
                         : mCategory == RawVariant::Type::INS ? "INS"sv
                         : mCategory == RawVariant::Type::DEL ? "DEL"sv
                         : mCategory == RawVariant::Type::MNP ? "MNP"sv
                                                              : "REF"sv;

  mInfoField = fmt::format(
      "{};{}TYPE={};LENGTH={};KMERLEN={}{}", vstate, is_str ? "STR;"sv : ""sv, vcategory, mVariantLength, kmerlen,
      is_str ? fmt::format(";STR_LEN={};STR_MOTIF={}", var->mStrResult.mStrLen, var->mStrResult.mStrMotif) : "");
}

auto VariantCall::AsVcfRecord() const -> std::string {
  // No newline. caller of this method will add new line if needed
  return fmt::format("{CHROM}\t{POS}\t.\t{REF}\t{ALT}\t{QUAL:.2f}\t.\t{INFO}\t{FORMAT}", fmt::arg("CHROM", mChromName),
                     fmt::arg("POS", mStartPos1), fmt::arg("REF", mRefAllele), fmt::arg("ALT", mAltAllele),
                     fmt::arg("QUAL", mSiteQuality), fmt::arg("INFO", mInfoField),
                     fmt::arg("FORMAT", absl::StrJoin(mFormatFields, "\t")));
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
