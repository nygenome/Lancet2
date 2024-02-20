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

  void AddEvidence(u32 rname_hash, Allele allele, Strand strand, u8 base_qual, u8 map_qual, u8 aln_diff_score);

  [[nodiscard]] inline auto RefFwdCount() const noexcept -> usize { return mRefFwdBaseQuals.size(); }
  [[nodiscard]] inline auto RefRevCount() const noexcept -> usize { return mRefRevBaseQuals.size(); }
  [[nodiscard]] inline auto AltFwdCount() const noexcept -> usize { return mAltFwdBaseQuals.size(); }
  [[nodiscard]] inline auto AltRevCount() const noexcept -> usize { return mAltRevBaseQuals.size(); }

  [[nodiscard]] inline auto TotalRefCov() const noexcept -> usize { return RefFwdCount() + RefRevCount(); }
  [[nodiscard]] inline auto TotalAltCov() const noexcept -> usize { return AltFwdCount() + AltRevCount(); }
  [[nodiscard]] inline auto TotalSampleCov() const noexcept -> usize { return TotalRefCov() + TotalAltCov(); }

  [[nodiscard]] auto AltFrequency() const -> f64;

  /// Normalized Phred scaled genotype likelihoods for all possible genotype combinations
  [[nodiscard]] auto ComputePLs() const -> std::array<int, 3>;

  /// Phred scaled quality median and range for the REF and ALT alleles
  struct QualStats {
    std::array<u8, 2> ref_alt_medians;
    std::array<u8, 2> ref_alt_ranges;
  };

  [[nodiscard]] auto AlleleQualityStats() const -> QualStats;
  [[nodiscard]] auto MappingQualityStats() const -> QualStats;
  [[nodiscard]] auto AlnDiffScoreStats() const -> QualStats;

 private:
  using Qualities = std::vector<u8>;
  using ReadNames = absl::flat_hash_map<u32, Strand>;

  ReadNames mRefNameHashes;
  ReadNames mAltNameHashes;

  Qualities mRefFwdBaseQuals;
  Qualities mRefRevBaseQuals;
  Qualities mAltFwdBaseQuals;
  Qualities mAltRevBaseQuals;

  Qualities mRefMapQuals;
  Qualities mAltMapQuals;

  Qualities mRefAlnDiffScores;
  Qualities mAltAlnDiffScores;

  [[nodiscard]] auto MeanErrorProbability(Allele allele) const -> f64;
  [[nodiscard]] auto BinomialSuccessRatios() const -> std::array<f64, 2>;
  [[nodiscard]] static auto ConvertGtProbsToPls(const std::array<f64, 3>& gt_probs) -> std::array<int, 3>;
  [[nodiscard]] static auto MedianOfSortedVector(absl::Span<const u8> data) -> u8;
  [[nodiscard]] static auto RangeOfSortedVector(absl::Span<const u8> data) -> u8;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
