#ifndef SRC_LANCET_BASE_COMPUTE_STATS_H_
#define SRC_LANCET_BASE_COMPUTE_STATS_H_

#include "lancet/base/types.h"

#include "absl/types/span.h"

#include <algorithm>
#include <concepts>
#include <limits>
#include <numeric>
#include <vector>

#include <cmath>
#include <cstdlib>

namespace lancet::base {

namespace detail {

template <std::floating_point T1 = f64, std::floating_point T2 = f64>
static constexpr auto AlmostEq(T1 first, T2 second) -> bool {
  auto const is_second_low_prec =
      std::numeric_limits<T1>::epsilon() > std::numeric_limits<T2>::epsilon();
  auto const a_val = is_second_low_prec ? first : static_cast<T2>(first);
  auto const b_val = is_second_low_prec ? static_cast<T1>(second) : second;
  auto const epsilon =
      is_second_low_prec ? std::numeric_limits<T1>::epsilon() : std::numeric_limits<T2>::epsilon();
  return std::abs(a_val - b_val) <= epsilon;
}

}  // namespace detail

template <typename T>
concept Number = std::integral<T> || std::floating_point<T>;

// ============================================================================
// OnlineStats — Welford's online algorithm for numerically stable statistics
//
// In plain terms: OnlineStats computes the mean and spread (standard
// deviation) of a stream of numbers — like base qualities or coverage
// values — without storing all the numbers. Each new number updates a
// running average and a running "how spread out are they?" counter.
// The clever part is that it avoids the naive formula (which can give
// wildly wrong answers when numbers are large and close together) by
// tracking deviations from the running mean instead. Two independent
// counters can be merged (e.g., from different threads) without
// re-reading the original data.
//
// WHY WELFORD'S: The naive formula Var = E[x²] − (E[x])² suffers from
// catastrophic cancellation when the mean is large relative to the spread
// (e.g., base qualities with mean ~35 and σ ~2).  Welford's recurrence
// maintains a running mean (m1) and sum-of-squared-deviations (m2)
// updated per observation, avoiding subtraction of nearly equal quantities.
//
// RECURRENCE (Add):
//   n ← n + 1
//   δ = x − m1                     // deviation from current mean
//   m1 ← m1 + δ/n                  // updated mean
//   m2 ← m2 + δ × (δ/n) × (n−1)   // Σ(xᵢ − mean)² accumulator
//
// FINAL STATISTICS:
//   Mean     = m1
//   Variance = m2 / (n − 1)    ← Bessel's correction (divides by n−1 instead
//                            of n to avoid underestimating spread from a sample)
//   StdDev   = √Variance
//
// PARALLEL MERGE (Chan et al.):
//   Two independent OnlineStats (n₁, m1₁, m2₁) and (n₂, m1₂, m2₂) merge as:
//   m2_combined = m2₁ + m2₂ + δ² × n₁ × n₂ / (n₁ + n₂)
//   where δ = m1₂ − m1₁
//
// PROPERTIES: O(1) memory, O(1) per update, single-pass, mergeable.
// ============================================================================
class OnlineStats {
 public:
  OnlineStats() = default;

  template <Number T>
  void Add(T const value) {
    auto const sample = static_cast<f64>(value);
    auto const old_num = mNum++;
    auto const delta = sample - mMoment1;
    auto const normalized_delta = delta / static_cast<f64>(mNum);

    mMoment1 += normalized_delta;
    mMoment2 += (delta * normalized_delta * static_cast<f64>(old_num));
  }

  void Clear() {
    mNum = 0;
    mMoment1 = 0.0;
    mMoment2 = 0.0;
  }

  void Merge(OnlineStats const& other) {
    auto const new_num = mNum + other.mNum;
    auto const delta1 = other.mMoment1 - mMoment1;
    auto const delta2 = delta1 * delta1;

    auto const fnum = static_cast<f64>(mNum);
    auto const other_fnum = static_cast<f64>(other.mNum);

    mMoment1 = ((fnum * mMoment1) + (other_fnum * other.mMoment1)) / static_cast<f64>(new_num);
    mMoment2 = mMoment2 + other.mMoment2 + (delta2 * fnum * other_fnum / static_cast<f64>(new_num));
    mNum = new_num;
  }

  [[nodiscard]] auto IsEmpty() const -> bool { return mNum == 0; }
  [[nodiscard]] auto Count() const -> usize { return mNum; }
  [[nodiscard]] auto Mean() const -> f64 { return mMoment1; }
  [[nodiscard]] auto Variance() const -> f64 {
    return mNum < 2 ? 0 : mMoment2 / static_cast<f64>(mNum - 1);
  }
  [[nodiscard]] auto StdDev() const -> f64 { return std::sqrt(Variance()); }

  auto operator==(OnlineStats const& rhs) const -> bool {
    return mNum == rhs.mNum &&
           detail::AlmostEq(mMoment1, rhs.mMoment1) &&
           detail::AlmostEq(mMoment2, rhs.mMoment2);
  }

  auto operator!=(OnlineStats const& rhs) const -> bool { return !(*this == rhs); }

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  usize mNum = 0;
  f64 mMoment1 = 0.0;
  f64 mMoment2 = 0.0;
};

template <Number T>
[[nodiscard]] inline auto Mean(absl::Span<T const> data) -> f64 {
  if (data.empty()) return 0.0;

  if (data.size() == 1) return data[0];

  static auto const SUMMER = [](f64 const sum, T const& num) -> f64 {
    return sum + static_cast<f64>(num);
  };
  f64 const sum = std::accumulate(data.cbegin(), data.cend(), 0.0, SUMMER);
  return sum / static_cast<f64>(data.size());
}

template <Number T>
[[nodiscard]] inline auto Median(absl::Span<T const> data) -> T {
  if (data.empty()) return 0;

  if (data.size() == 1) return data[0];

  std::vector<T> dcopy(data.cbegin(), data.cend());
  std::nth_element(dcopy.begin(), dcopy.begin() + (data.length() / 2), dcopy.end());
  T const half_item = dcopy[data.length() / 2];
  if (data.length() % 2 == 1) return half_item;

  std::nth_element(dcopy.begin(), dcopy.begin() + (data.length() / 2) - 1, dcopy.end());
  T const half_minus_one_item = dcopy[(data.length() / 2) - 1];
  return (half_item + half_minus_one_item) / 2;
}

template <Number T>
[[nodiscard]] inline auto Minimum(absl::Span<T const> data) -> T {
  if (data.empty()) return static_cast<T>(0);
  return *std::ranges::min_element(data);
}

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_COMPUTE_STATS_H_
