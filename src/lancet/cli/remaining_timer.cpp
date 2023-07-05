#include "lancet/cli/remaining_timer.h"

#include <cmath>

namespace lancet::cli {

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
RemainingTimer::RemainingTimer(const usize num_iterations, const usize num_threads)
    : mNumTotal(num_iterations), mNumThreads(num_threads) {
  static constexpr f64 PRIOR_MEAN_NS = 1e9;
  static constexpr f64 PRIOR_STDDEV_NS = 5e8;
  mRunStats.SetPriorMean(PRIOR_MEAN_NS);
  mRunStats.SetPriorStandardDeviation(PRIOR_STDDEV_NS);
}

void RemainingTimer::Update(const absl::Duration& loop_time) {
  mNumDone++;
  mRunStats.Add(absl::ToInt64Nanoseconds(loop_time));
}

auto RemainingTimer::EstimateRemaining() const -> absl::Duration {
  const auto estimated_ns_remaining = static_cast<f64>(mNumTotal - mNumDone) * mRunStats.Mean();
  return absl::Nanoseconds(estimated_ns_remaining / static_cast<f64>(mNumThreads));
}

auto RemainingTimer::MeanRatePerSecond() const -> f64 {
  static constexpr f64 NS_TO_SECS = 1e-9;
  static constexpr f64 WINDOWS_PER_SECOND_CONVERTER = -1.0;
  const auto ns_per_window = mRunStats.Mean() / static_cast<f64>(mNumThreads);
  return std::pow(ns_per_window * NS_TO_SECS, WINDOWS_PER_SECOND_CONVERTER);
}

}  // namespace lancet::cli
