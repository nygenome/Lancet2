#include "lancet/cli/remaining_timer.h"

namespace lancet::cli {

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
RemainingTimer::RemainingTimer(const usize num_iterations, const f64 prior_mean, const f64 prior_sd)
    : mNumTotal(num_iterations) {
  mRunStats.SetPriorMean(prior_mean);
  mRunStats.SetPriorStandardDeviation(prior_sd);
}

void RemainingTimer::Update(const absl::Duration& loop_time) {
  mNumDone++;
  mRunStats.Add(absl::ToInt64Nanoseconds(loop_time));
}

auto RemainingTimer::EstimateRemaining() const -> absl::Duration {
  const auto estimated_ns_remaining = static_cast<i64>(mNumTotal - mNumDone) * static_cast<i64>(mRunStats.Mean());
  return absl::Nanoseconds(estimated_ns_remaining);
}

}  // namespace lancet::cli
