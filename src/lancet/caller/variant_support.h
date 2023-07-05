#ifndef SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
#define SRC_LANCET_CALLER_VARIANT_SUPPORT_H_

#include <array>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"
#include "lancet/base/types.h"

namespace lancet::caller {

enum class Allele : bool { REF, ALT };
enum class Strand : bool { FWD, REV };

class VariantSupport {
 public:
  VariantSupport() = default;

  void AddEvidence(u32 rname_hash, Allele allele, Strand strand, u8 quality);

  [[nodiscard]] inline auto RefFwdCount() const noexcept -> usize { return mRefFwdQuals.size(); }
  [[nodiscard]] inline auto RefRevCount() const noexcept -> usize { return mRefRevQuals.size(); }
  [[nodiscard]] inline auto AltFwdCount() const noexcept -> usize { return mAltFwdQuals.size(); }
  [[nodiscard]] inline auto AltRevCount() const noexcept -> usize { return mAltRevQuals.size(); }

  [[nodiscard]] inline auto TotalRefCov() const noexcept -> usize { return RefFwdCount() + RefRevCount(); }
  [[nodiscard]] inline auto TotalAltCov() const noexcept -> usize { return AltFwdCount() + AltRevCount(); }
  [[nodiscard]] inline auto TotalSampleCov() const noexcept -> usize { return TotalRefCov() + TotalAltCov(); }

  [[nodiscard]] auto AltFrequency() const -> f64;

  /// Normalized Phred scaled genotype likelihoods for all possible genotype combinations
  [[nodiscard]] auto ComputePLs() const -> std::array<int, 3>;

  /// Phred scaled probability of the REF and ALT allele
  [[nodiscard]] auto MeanHaplotypeQualities() const -> std::array<u8, 2>;

 private:
  using Qualities = std::vector<u8>;
  using ReadNames = absl::flat_hash_map<u32, Strand>;

  ReadNames mRefNameHashes;
  ReadNames mAltNameHashes;
  Qualities mRefFwdQuals;
  Qualities mRefRevQuals;
  Qualities mAltFwdQuals;
  Qualities mAltRevQuals;

  [[nodiscard]] auto MeanErrorProbability(Allele allele) const -> f64;
  [[nodiscard]] auto BinomialSuccessRatios() const -> std::array<f64, 2>;
  [[nodiscard]] static auto ConvertGtProbsToPls(const std::array<f64, 3>& gt_probs) -> std::array<int, 3>;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
