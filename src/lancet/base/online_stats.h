#ifndef SRC_LANCET_BASE_ONLINE_STATS_H_
#define SRC_LANCET_BASE_ONLINE_STATS_H_

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdlib>
#include <limits>

#include "absl/types/span.h"
#include "lancet/base/types.h"

namespace detail {

template <std::floating_point T1 = f64, std::floating_point T2 = f64>
inline static constexpr auto AlmostEq(T1 first, T2 second) -> bool {
  const auto is_second_low_prec = std::numeric_limits<T1>::epsilon() > std::numeric_limits<T2>::epsilon();
  const auto aval = is_second_low_prec ? first : static_cast<T2>(first);
  const auto bval = is_second_low_prec ? static_cast<T1>(second) : second;
  const auto epsilon = is_second_low_prec ? std::numeric_limits<T1>::epsilon() : std::numeric_limits<T2>::epsilon();
  return std::abs(aval - bval) <= epsilon;
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

  template <Number T>
  [[nodiscard]] static auto Mean(absl::Span<const T> data) -> f64 {
    OnlineStats stats;
    std::ranges::for_each(data, [&stats](const T item) { stats.Add(item); });
    return stats.Mean();
  }

 private:
  usize mNum = 0;
  f64 mMoment1 = 0.0;
  f64 mMoment2 = 0.0;
};

#endif  // SRC_LANCET_BASE_ONLINE_STATS_H_