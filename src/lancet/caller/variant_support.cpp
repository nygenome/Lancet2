#include "lancet/caller/variant_support.h"

#include <algorithm>
#include <cmath>
#include <ranges>
#include <vector>

#include "boost/math/distributions/binomial.hpp"
#include "lancet/base/compute_stats.h"
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
  const auto alt_total = TotalAltCov();
  return alt_total == 0 ? 0.0 : static_cast<f64>(alt_total) / static_cast<f64>(TotalSampleCov());
}

auto VariantSupport::ComputePLs() const -> std::array<int, 3> {
  const auto nref = static_cast<f64>(TotalRefCov());
  const auto nalt = static_cast<f64>(TotalAltCov());
  const auto total = static_cast<f64>(TotalSampleCov());
  const auto fraction_ref = static_cast<f64>(nref) / static_cast<f64>(total);
  const auto fraction_alt = static_cast<f64>(nalt) / static_cast<f64>(total);

  // NOLINTBEGIN(readability-braces-around-statements)
  if (fraction_ref == 1.0 && fraction_alt == 0.0) return {0, hts::MAX_PHRED_SCORE, hts::MAX_PHRED_SCORE};
  if (fraction_ref == 0.0 && fraction_alt == 1.0) return {hts::MAX_PHRED_SCORE, hts::MAX_PHRED_SCORE, 0};
  // NOLINTEND(readability-braces-around-statements)

  f64 ref_hom_ll = 0.0;
  f64 het_alt_ll = 0.0;
  f64 alt_hom_ll = 0.0;

  const boost::math::binomial_distribution<f64> ref_dist(total, fraction_ref);
  const boost::math::binomial_distribution<f64> alt_dist(total, fraction_alt);

  // REF_HOM probability log likelihood
  const auto prob_all_refs = boost::math::pdf(ref_dist, total) * boost::math::pdf(alt_dist, 0);
  ref_hom_ll += prob_all_refs == 0.0 ? 0.0 : std::log10(prob_all_refs);

  // HET_ALT probability log likelihood
  const auto prob_nref_nalt = (boost::math::cdf(ref_dist, nref) - boost::math::cdf(ref_dist, 0)) *
                              (boost::math::cdf(alt_dist, nalt) - boost::math::cdf(alt_dist, 0));
  het_alt_ll += prob_nref_nalt == 0.0 ? 0.0 : std::log10(prob_nref_nalt);

  // ALT_HOM probability log likelihood
  const auto prob_all_alts = boost::math::pdf(ref_dist, 0) * boost::math::pdf(alt_dist, total);
  alt_hom_ll += prob_all_alts == 0.0 ? 0.0 : std::log10(prob_all_alts);

  // PL = -10 * log10(P(G|D)) where P(G|D) is the likelihood of genotype G given data D
  static constexpr f64 LL_TO_PHRED_MULTIPLIER = -10.0;
  ref_hom_ll = LL_TO_PHRED_MULTIPLIER * ref_hom_ll;
  het_alt_ll = LL_TO_PHRED_MULTIPLIER * het_alt_ll;
  alt_hom_ll = LL_TO_PHRED_MULTIPLIER * alt_hom_ll;

  // Subtract minimum PL value to get normalized PL values
  const f64 min_ll_value = std::min({ref_hom_ll, het_alt_ll, alt_hom_ll});
  ref_hom_ll -= min_ll_value;
  het_alt_ll -= min_ll_value;
  alt_hom_ll -= min_ll_value;

  return {hts::ClampPhredScore(ref_hom_ll), hts::ClampPhredScore(het_alt_ll), hts::ClampPhredScore(alt_hom_ll)};
}

auto VariantSupport::StrandBiasScore() const -> u8 {
  using Row = hts::FisherExact::Row;
  const auto fwds = Row{static_cast<int>(mRefFwdQuals.size()), static_cast<int>(mAltFwdQuals.size())};
  const auto revs = Row{static_cast<int>(mRefRevQuals.size()), static_cast<int>(mAltRevQuals.size())};
  const auto result = hts::FisherExact::Test({fwds, revs});
  return hts::ErrorProbToPhred(result.mDiffProb);
}

}  // namespace lancet::caller
