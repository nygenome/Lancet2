#ifndef SRC_LANCET_BASE_COMPUTE_STATS_H_
#define SRC_LANCET_BASE_COMPUTE_STATS_H_

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdlib>
#include <limits>
#include <numeric>

#include "absl/types/span.h"
#include "lancet/base/types.h"

namespace detail {

template <std::floating_point T1 = f64, std::floating_point T2 = f64>
inline static constexpr auto AlmostEq(T1 first, T2 second) -> bool {
  const auto IsSecondLowPrec = std::numeric_limits<T1>::epsilon() > std::numeric_limits<T2>::epsilon();
  const auto Aval = IsSecondLowPrec ? first : static_cast<T2>(first);
  const auto Bval = IsSecondLowPrec ? static_cast<T1>(second) : second;
  const auto Epsilon = IsSecondLowPrec ? std::numeric_limits<T1>::epsilon() : std::numeric_limits<T2>::epsilon();
  return std::abs(Aval - Bval) <= Epsilon;
}

}  // namespace detail

template <typename T>
concept Number = std::integral<T> || std::floating_point<T>;

class OnlineStats {
 public:
  OnlineStats() = default;

  template <Number T>
  void Add(const T value) {
    const auto sample = static_cast<f64>(value);
    const auto old_num = mNum++;
    const auto delta = sample - mMoment1;
    const auto normalized_delta = delta / static_cast<f64>(mNum);

    mMoment1 += normalized_delta;
    mMoment2 += (delta * normalized_delta * static_cast<f64>(old_num));
  }

  void Clear() {
    mNum = 0;
    mMoment1 = 0.0;
    mMoment2 = 0.0;
  }

  void Merge(const OnlineStats& other) {
    const auto new_num = mNum + other.mNum;
    const auto delta1 = other.mMoment1 - mMoment1;
    const auto delta2 = delta1 * delta1;

    const auto fnum = static_cast<f64>(mNum);
    const auto other_fnum = static_cast<f64>(other.mNum);

    mMoment1 = ((fnum * mMoment1) + (other_fnum * other.mMoment1)) / static_cast<f64>(new_num);
    mMoment2 = mMoment2 + other.mMoment2 + delta2 * fnum * other_fnum / static_cast<f64>(new_num);
    mNum = new_num;
  }

  [[nodiscard]] auto IsEmpty() const -> bool { return mNum == 0; }
  [[nodiscard]] auto Count() const -> usize { return mNum; }
  [[nodiscard]] auto Mean() const -> f64 { return mMoment1; }
  [[nodiscard]] auto Variance() const -> f64 { return mNum < 2 ? 0 : mMoment2 / static_cast<f64>(mNum - 1); }
  [[nodiscard]] auto StdDev() const -> f64 { return std::sqrt(Variance()); }

  auto operator==(const OnlineStats& rhs) const -> bool {
    return mNum == rhs.mNum && detail::AlmostEq(mMoment1, rhs.mMoment1) && detail::AlmostEq(mMoment2, rhs.mMoment2);
  }

  auto operator!=(const OnlineStats& rhs) const -> bool { return !(*this == rhs); }

 private:
  usize mNum = 0;
  f64 mMoment1 = 0.0;
  f64 mMoment2 = 0.0;
};

template <Number T>
[[nodiscard]] inline auto Mean(absl::Span<const T> data) -> f64 {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (data.empty()) return 0.0;

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (data.size() == 1) return data[0];

  static const auto summer = [](const f64 sum, const T& num) { return sum + static_cast<f64>(num); };
  const f64 sum = std::accumulate(data.cbegin(), data.cend(), 0.0, summer);
  return sum / static_cast<f64>(data.size());
}

template <Number T>
[[nodiscard]] inline auto Median(absl::Span<const T> data) -> T {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (data.empty()) return 0;

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (data.size() == 1) return data[0];

  std::vector<T> dcopy(data.cbegin(), data.cend());
  std::nth_element(dcopy.begin(), dcopy.begin() + data.length() / 2, dcopy.end());
  const T half_item = dcopy[data.length() / 2];
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (data.length() % 2 == 1) return half_item;

  std::nth_element(dcopy.begin(), dcopy.begin() + (data.length() / 2) - 1, dcopy.end());
  const T half_minus_one_item = dcopy[(data.length() / 2) - 1];
  return (half_item + half_minus_one_item) / 2;
}

template <Number T>
[[nodiscard]] inline auto Minimum(absl::Span<const T> data) -> T {
  static const auto accumulator = [](const T curr_min, const T value) -> T { return std::min(curr_min, value); };
  const auto result = std::accumulate(data.cbegin(), data.cend(), std::numeric_limits<T>::max(), accumulator);
  return result == std::numeric_limits<T>::max() ? static_cast<T>(0) : result;
}

#endif  // SRC_LANCET_BASE_COMPUTE_STATS_H_
