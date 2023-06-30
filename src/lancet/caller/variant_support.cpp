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
  const auto nref = static_cast<f64>(TotalRefCov());
  const auto nalt = static_cast<f64>(TotalAltCov());
  const auto total = static_cast<f64>(TotalSampleCov());
  // NOLINTBEGIN(readability-braces-around-statements)
  if (nref == 0.0 && nalt == 0.0 && total == 0.0) return {0, 0, 0};
  if (nref == 0.0) return {hts::MAX_PHRED_SCORE, hts::MAX_PHRED_SCORE, 0};
  if (nalt == 0.0) return {0, hts::MAX_PHRED_SCORE, hts::MAX_PHRED_SCORE};
  // NOLINTEND(readability-braces-around-statements)

  const auto ref_fraction = static_cast<f64>(nref) / static_cast<f64>(total);
  const auto alt_fraction = static_cast<f64>(nalt) / static_cast<f64>(total);
  const boost::math::binomial_distribution<f64> ref_dist(total, ref_fraction);
  const boost::math::binomial_distribution<f64> alt_dist(total, alt_fraction);

  const auto prob_hom_ref = std::max(boost::math::pdf(ref_dist, total), boost::math::pdf(alt_dist, 0));
  const auto prob_het_alt = std::max(boost::math::pdf(ref_dist, nref), boost::math::pdf(alt_dist, nalt));
  const auto prob_hom_alt = std::max(boost::math::pdf(ref_dist, 0), boost::math::pdf(alt_dist, total));

  return ConvertGtProbsToPls({prob_hom_ref, prob_het_alt, prob_hom_alt});
}

auto VariantSupport::StrandBiasScore() const -> u8 {
  using Row = hts::FisherExact::Row;
  const auto fwds = Row{static_cast<int>(mRefFwdQuals.size()), static_cast<int>(mAltFwdQuals.size())};
  const auto revs = Row{static_cast<int>(mRefRevQuals.size()), static_cast<int>(mAltRevQuals.size())};
  const auto result = hts::FisherExact::Test({fwds, revs});
  return hts::ErrorProbToPhred(result.mDiffProb);
}

auto VariantSupport::MeanHaplotypeQualities() const -> std::array<int, 2> {
  const std::array<absl::Span<const u8>, 2> ref_quals{mRefFwdQuals, mRefRevQuals};
  const std::array<absl::Span<const u8>, 2> alt_quals{mAltFwdQuals, mAltRevQuals};
  const auto ref_range = std::ranges::join_view(ref_quals);
  const auto alt_range = std::ranges::join_view(alt_quals);

  static const auto summer = [](const f64 sum, const u8 bqual) -> f64 { return sum + hts::PhredToErrorProb(bqual); };
  const auto sum_ref_eprob = std::accumulate(ref_range.begin(), ref_range.end(), 0.0, summer);
  const auto sum_alt_eprob = std::accumulate(alt_range.begin(), alt_range.end(), 0.0, summer);
  const auto mean_ref_eprob = sum_ref_eprob == 0.0 ? 0.0 : sum_ref_eprob / static_cast<f64>(TotalRefCov());
  const auto mean_alt_eprob = sum_alt_eprob == 0.0 ? 0.0 : sum_alt_eprob / static_cast<f64>(TotalAltCov());

  return {hts::ErrorProbToPhred(mean_ref_eprob), hts::ErrorProbToPhred(mean_alt_eprob)};
}

auto VariantSupport::ConvertGtProbsToPls(const std::array<f64, 3>& gt_probs) -> std::array<int, 3> {
  const auto [prob_hom_ref, prob_het_alt, prob_hom_alt] = gt_probs;

  const f64 hom_ref_ll = prob_hom_ref == 0.0 ? 0.0 : std::log10(prob_hom_ref);
  const f64 het_alt_ll = prob_het_alt == 0.0 ? 0.0 : std::log10(prob_het_alt);
  const f64 hom_alt_ll = prob_hom_alt == 0.0 ? 0.0 : std::log10(prob_hom_alt);

  static constexpr f64 LL_TO_PHRED_MULTIPLIER = -10.0;
  auto hom_ref_phred = LL_TO_PHRED_MULTIPLIER * hom_ref_ll;
  auto het_alt_phred = LL_TO_PHRED_MULTIPLIER * het_alt_ll;
  auto hom_alt_phred = LL_TO_PHRED_MULTIPLIER * hom_alt_ll;

  // Subtract minimum phred value to get normalized PL values
  const f64 min_phred = std::min({hom_ref_phred, het_alt_phred, hom_alt_phred});
  const auto pl_hom_ref = hom_ref_phred - min_phred;
  const auto pl_het_alt = het_alt_phred - min_phred;
  const auto pl_hom_alt = hom_alt_phred - min_phred;

  return {hts::ClampPhredScore(pl_hom_ref), hts::ClampPhredScore(pl_het_alt), hts::ClampPhredScore(pl_hom_alt)};
}

}  // namespace lancet::caller
