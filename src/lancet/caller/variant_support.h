#ifndef SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
#define SRC_LANCET_CALLER_VARIANT_SUPPORT_H_

#include <array>
#include <cmath>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"
#include "boost/math/statistics/univariate_statistics.hpp"
#include "lancet/base/types.h"

namespace lancet::caller {

enum class Allele : bool { REF, ALT };
enum class Strand : bool { FWD, REV };

template <typename T>
concept Number = std::integral<T> || std::floating_point<T>;

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

  /// Minimum, Median, Maximum and Median Absolute Deviation for the REF and ALT alleles
  struct Statistics {
    int refMinVal;
    int refMedian;
    int refMaxVal;
    int refMADVal;

    int altMinVal;
    int altMedian;
    int altMaxVal;
    int altMADVal;
  };

  [[nodiscard]] auto AlleleQualityStats() const -> Statistics;
  [[nodiscard]] auto MappingQualityStats() const -> Statistics;
  [[nodiscard]] auto AlnDiffScoreStats() const -> Statistics;

 private:
  using Qualities = std::vector<u8>;
  using ReadNames = absl::flat_hash_map<u32, Strand>;

  ReadNames mRefNameHashes{};
  ReadNames mAltNameHashes{};

  Qualities mRefFwdBaseQuals{};
  Qualities mRefRevBaseQuals{};
  Qualities mAltFwdBaseQuals{};
  Qualities mAltRevBaseQuals{};

  Qualities mRefMapQuals{};
  Qualities mAltMapQuals{};

  Qualities mRefAlnDiffScores{};
  Qualities mAltAlnDiffScores{};

  [[nodiscard]] auto MeanErrorProbability(Allele allele) const -> f64;
  [[nodiscard]] auto BinomialSuccessRatios() const -> std::array<f64, 2>;
  [[nodiscard]] static auto ConvertGtProbsToPls(const std::array<f64, 3>& gt_probs) -> std::array<int, 3>;

  template <Number T>
  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  [[nodiscard]] static auto BuildStats(absl::Span<const T> data_ref, absl::Span<const T> data_alt) -> Statistics {
    std::vector<f64> refs;
    std::vector<f64> alts;
    refs.reserve(data_ref.size());
    alts.reserve(data_alt.size());

    std::ranges::for_each(data_ref, [&refs](const u8 bqual) { refs.push_back(static_cast<f64>(bqual)); });
    std::ranges::for_each(data_alt, [&alts](const u8 bqual) { alts.push_back(static_cast<f64>(bqual)); });
    std::ranges::sort(refs);
    std::ranges::sort(alts);

    const auto sz_ref = refs.size();
    const auto sz_alt = alts.size();

    const auto min_ref = refs.empty() ? 0 : static_cast<int>(refs[0]);
    const auto min_alt = alts.empty() ? 0 : static_cast<int>(alts[0]);
    const auto max_ref = refs.empty() ? 0 : static_cast<int>(refs[sz_ref - 1]);
    const auto max_alt = alts.empty() ? 0 : static_cast<int>(alts[sz_alt - 1]);

    const auto ref_median = refs.empty()        ? 0.0
                            : (sz_ref % 2 == 0) ? (refs[sz_ref / 2 - 1] + refs[sz_ref / 2]) / 2.0
                                                : refs[sz_ref / 2];

    const auto alt_median = alts.empty()        ? 0.0
                            : (sz_alt % 2 == 0) ? (alts[sz_alt / 2 - 1] + alts[sz_alt / 2]) / 2.0
                                                : alts[sz_alt / 2];

    const auto ref_mad = refs.empty() ? 0.0 : boost::math::statistics::median_absolute_deviation(refs, ref_median);
    const auto alt_mad = alts.empty() ? 0.0 : boost::math::statistics::median_absolute_deviation(alts, alt_median);

    return {
        .refMinVal = min_ref,
        .refMedian = static_cast<int>(std::round(ref_median)),
        .refMaxVal = max_ref,
        .refMADVal = static_cast<int>(std::round(ref_mad)),

        .altMinVal = min_alt,
        .altMedian = static_cast<int>(std::round(alt_median)),
        .altMaxVal = max_alt,
        .altMADVal = static_cast<int>(std::round(alt_mad)),
    };
  }
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
