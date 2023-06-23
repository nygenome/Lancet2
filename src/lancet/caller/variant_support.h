#ifndef SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
#define SRC_LANCET_CALLER_VARIANT_SUPPORT_H_

#include <array>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/types/span.h"
#include "lancet/base/types.h"

namespace lancet::caller {

enum class Allele : bool { REF, ALT };
enum class Strand : bool { FWD, REV };
using AlleleStrand = std::pair<Allele, Strand>;

class VariantSupport {
 public:
  VariantSupport() = default;

  void AddQual(u8 qual, const AlleleStrand& allele_strand);

  [[nodiscard]] inline auto RefFwdCount() const noexcept -> usize { return mRefFwdQuals.size(); }
  [[nodiscard]] inline auto RefRevCount() const noexcept -> usize { return mRefRevQuals.size(); }
  [[nodiscard]] inline auto AltFwdCount() const noexcept -> usize { return mAltFwdQuals.size(); }
  [[nodiscard]] inline auto AltRevCount() const noexcept -> usize { return mAltRevQuals.size(); }

  [[nodiscard]] inline auto TotalRefCov() const noexcept -> usize { return RefFwdCount() + RefRevCount(); }
  [[nodiscard]] inline auto TotalAltCov() const noexcept -> usize { return AltFwdCount() + AltRevCount(); }
  [[nodiscard]] inline auto TotalSampleCov() const noexcept -> usize { return TotalRefCov() + TotalAltCov(); }

  [[nodiscard]] auto AltFrequency() const -> f64;

  /// Normalized Phred scaled genotype likelihoods for all possible genotype combinations
  [[nodiscard]] auto NormalizedPhredLikelihoods() const -> std::array<int, 3>;

  /// Phred scaled probability of strand bias being present in the ref and alt alleles
  [[nodiscard]] auto StrandBiasScore() const -> u8;

 private:
  using QualsVector = std::vector<u8>;
  QualsVector mRefFwdQuals;
  QualsVector mRefRevQuals;
  QualsVector mAltFwdQuals;
  QualsVector mAltRevQuals;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
