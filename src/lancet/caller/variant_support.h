#ifndef SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
#define SRC_LANCET_CALLER_VARIANT_SUPPORT_H_

#include <array>
#include <concepts>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "lancet/base/types.h"

namespace lancet::caller {

enum class Allele : bool { REF, ALT };
enum class Strand : bool { FWD, REV };

template <typename T>
concept Number = std::integral<T> || std::floating_point<T>;

class VariantSupport {
 public:
  VariantSupport() = default;

  void AddEvidence(u32 rname_hash, Allele allele, Strand strand, u8 base_qual);

  [[nodiscard]] auto RefFwdCount() const noexcept -> usize { return mRefFwdBaseQuals.size(); }
  [[nodiscard]] auto RefRevCount() const noexcept -> usize { return mRefRevBaseQuals.size(); }
  [[nodiscard]] auto AltFwdCount() const noexcept -> usize { return mAltFwdBaseQuals.size(); }
  [[nodiscard]] auto AltRevCount() const noexcept -> usize { return mAltRevBaseQuals.size(); }

  [[nodiscard]] auto TotalRefCov() const noexcept -> usize { return RefFwdCount() + RefRevCount(); }
  [[nodiscard]] auto TotalAltCov() const noexcept -> usize { return AltFwdCount() + AltRevCount(); }
  [[nodiscard]] auto TotalSampleCov() const noexcept -> usize { return TotalRefCov() + TotalAltCov(); }


  /// Normalized Phred scaled genotype likelihoods for all possible genotype combinations
  [[nodiscard]] auto ComputePLs() const -> std::array<int, 3>;

 private:
  using Qualities = std::vector<u8>;
  using ReadNames = absl::flat_hash_map<u32, Strand>;

  ReadNames mRefNameHashes;
  ReadNames mAltNameHashes;

  Qualities mRefFwdBaseQuals;
  Qualities mRefRevBaseQuals;
  Qualities mAltFwdBaseQuals;
  Qualities mAltRevBaseQuals;

  [[nodiscard]] auto MeanErrorProbability(Allele allele) const -> f64;
  [[nodiscard]] auto BinomialSuccessRatios() const -> std::array<f64, 2>;
  [[nodiscard]] static auto ConvertGtProbsToPls(const std::array<f64, 3>& gt_probs) -> std::array<int, 3>;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
