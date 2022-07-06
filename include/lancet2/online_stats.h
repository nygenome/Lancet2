#pragma once

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>

namespace lancet2 {
namespace detail {
template <typename T = float>
inline static constexpr auto AlmostEqual(T a, T b) -> bool {
  return std::abs(a - b) <= std::numeric_limits<T>::epsilon();
}

template <typename T1 = float, typename T2 = float>
inline static constexpr auto AlmostEqual(T1 a, T2 b) -> bool {
  const auto isALowPrec = std::numeric_limits<T1>::epsilon() > std::numeric_limits<T2>::epsilon();
  const auto aval = isALowPrec ? a : static_cast<T2>(a);
  const auto bval = isALowPrec ? static_cast<T1>(b) : b;
  const auto epsilon = isALowPrec ? std::numeric_limits<T1>::epsilon() : std::numeric_limits<T2>::epsilon();
  return std::abs(aval - bval) <= epsilon;
}
}  // namespace detail

class OnlineStats {
 public:
  OnlineStats() = default;

  void Add(const double sample) {
    const auto oldNum = num++;
    const auto delta = sample - moment1;
    const auto normalizedDelta = delta / static_cast<double>(num);

    moment1 += normalizedDelta;
    moment2 += (delta * normalizedDelta * static_cast<double>(oldNum));
  }

  void Clear() {
    num = 0;
    moment1 = 0.0;
    moment2 = 0.0;
  }

  auto Merge(const OnlineStats& other) -> OnlineStats& {
    const auto newNum = num + other.num;
    const auto delta1 = other.moment1 - moment1;
    const auto delta2 = delta1 * delta1;

    moment1 = (( static_cast<double>(num) * moment1) + (static_cast<double>(other.num) * other.moment1)) / static_cast<double>(newNum);
    moment2 = moment2 + other.moment2 + delta2 * static_cast<double>(num) * static_cast<double>(other.num) / static_cast<double>(newNum);
    num = newNum;

    return *this;
  }

  [[nodiscard]] auto IsEmpty() const -> bool { return num == 0; }
  [[nodiscard]] auto Size() const -> std::size_t { return num; }
  [[nodiscard]] auto Mean() const -> double { return moment1; }
  [[nodiscard]] auto Variance() const -> double { return num < 2 ? 0 : moment2 / static_cast<double>(num - 1); }
  [[nodiscard]] auto StandardDeviation() const -> double { return std::sqrt(Variance()); }

  auto operator==(const OnlineStats& rhs) const -> bool {
    return num == rhs.num && detail::AlmostEqual(moment1, rhs.moment1) && detail::AlmostEqual(moment2, rhs.moment2);
  }

  auto operator!=(const OnlineStats& rhs) const -> bool { return !(*this == rhs); }

  friend auto operator+(const OnlineStats& lhs, const OnlineStats& rhs) -> OnlineStats {
    OnlineStats result;
    result = result.Merge(lhs);
    result = result.Merge(rhs);
    return result;
  }

 private:
  std::size_t num = 0;
  double moment1 = 0.0;
  double moment2 = 0.0;
};
}  // namespace lancet2
