#include "lancet/caller/variant_support.h"

#include <algorithm>
#include <cmath>
#include <ranges>

#include "absl/container/fixed_array.h"
#include "lancet/base/compute_stats.h"
#include "lancet/hts/fisher_exact.h"

namespace lancet::caller {

void VariantSupport::AddEvidence(const u32 rname_hash, const Allele allele, const Strand strand, const u8 quality) {
  switch (allele) {
    case Allele::REF:
      mRefNameHashes.emplace(rname_hash);
      strand == Strand::FWD ? mRefFwdQuals.emplace_back(quality) : mRefRevQuals.emplace_back(quality);
      break;
    default:
      mAltNameHashes.emplace(rname_hash);
      strand == Strand::FWD ? mAltFwdQuals.emplace_back(quality) : mAltRevQuals.emplace_back(quality);
      break;
  }
}

auto VariantSupport::AltFrequency() const -> f64 {
  const auto alt_total = TotalAltCov();
  return alt_total == 0 ? 0.0 : static_cast<f64>(alt_total) / static_cast<f64>(TotalSampleCov());
}

auto VariantSupport::NormalizedPhredLikelihoods() const -> std::array<int, 3> {
  const std::array<absl::Span<const u8>, 2> ref_quals{mRefFwdQuals, mRefRevQuals};
  const std::array<absl::Span<const u8>, 2> alt_quals{mAltFwdQuals, mAltRevQuals};

  std::vector<f64> ref_probs;
  std::vector<f64> alt_probs;
  ref_probs.reserve(TotalRefCov());
  alt_probs.reserve(TotalAltCov());

  std::ranges::transform(std::ranges::join_view(ref_quals), std::back_inserter(ref_probs), hts::PhredToErrorProb);
  std::ranges::transform(std::ranges::join_view(alt_quals), std::back_inserter(alt_probs), hts::PhredToErrorProb);

  // Assume bi-allelic diploid calls
  static constexpr f64 HET_RATIO = 0.5;
  f64 ref_hom_ll = 0.0;
  f64 het_alt_ll = 0.0;
  f64 alt_hom_ll = 0.0;

  for (const f64 ref_ep : ref_probs) {
    ref_hom_ll += ref_ep == 1.0 ? 0.0 : std::log10(1.0 - ref_ep);               // RR
    het_alt_ll += std::log10(HET_RATIO * (1.0 - ref_ep) + HET_RATIO * ref_ep);  // RA
    alt_hom_ll += std::log10(ref_ep);                                           // AA
  }

  for (const f64 alt_ep : alt_probs) {
    ref_hom_ll += std::log10(alt_ep);                                           // RR
    het_alt_ll += std::log10(HET_RATIO * (1.0 - alt_ep) + HET_RATIO * alt_ep);  // AR
    alt_hom_ll += alt_ep == 1.0 ? 0.0 : std::log10(1.0 - alt_ep);               // AA
  }

  // PL = -10 * log10(P(G|D)) where P(G|D) is the likelihood of genotype G given data D
  static constexpr f64 LL_TO_PHRED_MULTIPLIER = -10.0;
  ref_hom_ll = LL_TO_PHRED_MULTIPLIER * ref_hom_ll;
  het_alt_ll = LL_TO_PHRED_MULTIPLIER * het_alt_ll;
  alt_hom_ll = LL_TO_PHRED_MULTIPLIER * alt_hom_ll;

  // Subtract minimum PL value to get normalized PL values
  const f64 min_ll_value = std::min({ref_hom_ll, het_alt_ll, alt_hom_ll});
  return {static_cast<int>(std::round(ref_hom_ll - min_ll_value)),
          static_cast<int>(std::round(het_alt_ll - min_ll_value)),
          static_cast<int>(std::round(alt_hom_ll - min_ll_value))};
}

auto VariantSupport::StrandBiasScore() const -> u8 {
  using Row = hts::FisherExact::Row;
  const auto fwds = Row{static_cast<int>(mRefFwdQuals.size()), static_cast<int>(mAltFwdQuals.size())};
  const auto revs = Row{static_cast<int>(mRefRevQuals.size()), static_cast<int>(mAltRevQuals.size())};
  const auto result = hts::FisherExact::Test({fwds, revs});
  return hts::ErrorProbToPhred(result.mDiffProb);
}

auto VariantSupport::SupportingReadHashes(const Allele allele) const -> roaring::Roaring {
  const auto& curr_allele_set = allele == Allele::REF ? mRefNameHashes : mAltNameHashes;
  const absl::FixedArray<u32> allele_hashes(curr_allele_set.cbegin(), curr_allele_set.cend());
  return {curr_allele_set.size(), allele_hashes.data()};
}

}  // namespace lancet::caller
