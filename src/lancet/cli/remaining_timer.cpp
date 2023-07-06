#include "lancet/cli/remaining_timer.h"

#include <cmath>

namespace lancet::cli {

RemainingTimer::RemainingTimer(const usize num_iterations) : mNumTotal(num_iterations) {}

void RemainingTimer::Increment() {
  mNumDone++;
  mRunStats.Add(absl::ToInt64Nanoseconds(mProgressTimer.Runtime()));
}

auto RemainingTimer::EstimateRemaining() const -> absl::Duration {
  const auto estimated_ns_remaining = static_cast<f64>(mNumTotal - mNumDone) * mRunStats.Mean();
  return absl::Nanoseconds(estimated_ns_remaining);
}

auto RemainingTimer::MeanRatePerSecond() const -> f64 {
  static constexpr f64 NS_TO_SECS = 1e-9;
  static constexpr f64 WINDOWS_PER_SECOND_CONVERTER = -1.0;
  return std::pow(mRunStats.Mean() * NS_TO_SECS, WINDOWS_PER_SECOND_CONVERTER);
}

}  // namespace lancet::cli
