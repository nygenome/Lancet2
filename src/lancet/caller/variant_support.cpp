#include "lancet/caller/variant_support.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ranges>

#include "boost/math/distributions/binomial.hpp"
#include "lancet/hts/fisher_exact.h"

namespace lancet::caller {

void VariantSupport::AddEvidence(const u32 rname_hash, const Allele allele, const Strand strand, const u8 quality) {
  const auto& name_hashes = allele == Allele::REF ? mRefNameHashes : mAltNameHashes;
  const auto itr = name_hashes.find(rname_hash);
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (itr != name_hashes.end() && itr->second == strand) return;

  switch (allele) {
    case Allele::REF:
      mRefNameHashes.emplace(rname_hash, strand);
      strand == Strand::FWD ? mRefFwdQuals.emplace_back(quality) : mRefRevQuals.emplace_back(quality);
      break;
    default:
      mAltNameHashes.emplace(rname_hash, strand);
      strand == Strand::FWD ? mAltFwdQuals.emplace_back(quality) : mAltRevQuals.emplace_back(quality);
      break;
  }
}

auto VariantSupport::AltFrequency() const -> f64 {
  return TotalAltCov() == 0 ? 0.0 : static_cast<f64>(TotalAltCov()) / static_cast<f64>(TotalSampleCov());
}

auto VariantSupport::ComputePLs() const -> std::array<int, 3> {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (TotalSampleCov() == 0) return {0, 0, 0};

  const auto total_count = static_cast<f64>(TotalSampleCov());
  const auto [success_ratio_ref, success_ratio_alt] = BinomialSuccessRatios();
  const boost::math::binomial_distribution<f64> ref_dist(total_count, success_ratio_ref);
  const boost::math::binomial_distribution<f64> alt_dist(total_count, success_ratio_alt);
  const auto prob_hom_ref = boost::math::pdf(ref_dist, total_count);
  const auto prob_hom_alt = boost::math::pdf(alt_dist, total_count);
  const auto prob_het_alt = 1.0 - (prob_hom_ref + prob_hom_alt);

  return ConvertGtProbsToPls({prob_hom_ref, prob_het_alt, prob_hom_alt});
}

auto VariantSupport::StrandBiasScore() const -> u8 {
  using Row = hts::FisherExact::Row;
  const auto fwds = Row{static_cast<int>(mRefFwdQuals.size()), static_cast<int>(mAltFwdQuals.size())};
  const auto revs = Row{static_cast<int>(mRefRevQuals.size()), static_cast<int>(mAltRevQuals.size())};
  const auto result = hts::FisherExact::Test({fwds, revs});
  return hts::ErrorProbToPhred(result.mDiffProb);
}

auto VariantSupport::MeanHaplotypeQualities() const -> std::array<u8, 2> {
  return {hts::ErrorProbToPhred(MeanErrorProbability(Allele::REF)),
          hts::ErrorProbToPhred(MeanErrorProbability(Allele::ALT))};
}

auto VariantSupport::MeanErrorProbability(const Allele allele) const -> f64 {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (TotalSampleCov() == 0) return std::numeric_limits<f32>::min();

  const auto total_allele_cov = allele == Allele::REF ? TotalRefCov() : TotalAltCov();
  const auto data = allele == Allele::REF ? std::array<absl::Span<const u8>, 2>{mRefFwdQuals, mRefRevQuals}
                                          : std::array<absl::Span<const u8>, 2>{mAltFwdQuals, mAltRevQuals};

  const auto quals = std::ranges::join_view(data);
  static const auto summer = [](const f64 sum, const u8 bql) { return sum + hts::PhredToErrorProb(bql); };
  const auto err_prob_sum = std::accumulate(quals.begin(), quals.end(), 0.0, summer);
  return err_prob_sum == 0.0 ? std::numeric_limits<f32>::min() : err_prob_sum / static_cast<f64>(total_allele_cov);
}

auto VariantSupport::BinomialSuccessRatios() const -> std::array<f64, 2> {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (TotalSampleCov() == 0) return {std::numeric_limits<f32>::min(), std::numeric_limits<f32>::min()};

  const auto ref_count = static_cast<f64>(TotalRefCov());
  const auto alt_count = static_cast<f64>(TotalAltCov());
  const auto total_count = static_cast<f64>(TotalSampleCov());
  const auto ref_err_prob = MeanErrorProbability(Allele::REF);
  const auto alt_err_prob = MeanErrorProbability(Allele::ALT);

  static constexpr f64 MIN_PICK_PROB = 0.0 + std::numeric_limits<f32>::min();
  static constexpr f64 MAX_PICK_PROB = 1.0 - std::numeric_limits<f32>::min();

  if (alt_count == 0.0) {
    // REF allele is the most likely allele to be picked if alt_count == 0.0
    const auto prob_pick_ref = std::clamp(1.0 - (ref_err_prob / total_count), MIN_PICK_PROB, MAX_PICK_PROB);
    const auto prob_pick_alt = 1.0 - prob_pick_ref;
    return {prob_pick_ref, prob_pick_alt};
  }

  if (ref_count == 0.0) {
    // ALT allele is the most likely allele to be picked if ref_count == 0.0
    const auto prob_pick_alt = std::clamp(1.0 - (alt_err_prob / total_count), MIN_PICK_PROB, MAX_PICK_PROB);
    const auto prob_pick_ref = 1.0 - prob_pick_alt;
    return {prob_pick_ref, prob_pick_alt};
  }

  const auto weight_ref = std::clamp((1.0 - ref_err_prob) + alt_err_prob, 0.0, 1.0);
  const auto weight_alt = std::clamp((1.0 - alt_err_prob) + ref_err_prob, 0.0, 1.0);
  const auto prob_pick_ref = (ref_count / total_count) * weight_ref;
  const auto prob_pick_alt = (alt_count / total_count) * weight_alt;
  return {prob_pick_ref, prob_pick_alt};
}

auto VariantSupport::ConvertGtProbsToPls(const std::array<f64, 3>& gt_probs) -> std::array<int, 3> {
  const auto [prob_hom_ref, prob_het_alt, prob_hom_alt] = gt_probs;

  static constexpr f64 min_log_prob = std::numeric_limits<f64>::min_exponent10;
  const f64 hom_ref_ll = prob_hom_ref == 0.0 ? min_log_prob : std::log10(prob_hom_ref);
  const f64 het_alt_ll = prob_het_alt == 0.0 ? min_log_prob : std::log10(prob_het_alt);
  const f64 hom_alt_ll = prob_hom_alt == 0.0 ? min_log_prob : std::log10(prob_hom_alt);

  static constexpr f64 LL_TO_PHRED_MULTIPLIER = -10.0;
  auto hom_ref_phred = LL_TO_PHRED_MULTIPLIER * hom_ref_ll;
  auto het_alt_phred = LL_TO_PHRED_MULTIPLIER * het_alt_ll;
  auto hom_alt_phred = LL_TO_PHRED_MULTIPLIER * hom_alt_ll;

  // Subtract minimum phred value to get normalized PL values
  const f64 min_phred = std::min({hom_ref_phred, het_alt_phred, hom_alt_phred});
  const auto pl_hom_ref = std::ceil(hom_ref_phred - min_phred);
  const auto pl_het_alt = std::ceil(het_alt_phred - min_phred);
  const auto pl_hom_alt = std::ceil(hom_alt_phred - min_phred);

  return {static_cast<int>(pl_hom_ref), static_cast<int>(pl_het_alt), static_cast<int>(pl_hom_alt)};
}

}  // namespace lancet::caller
