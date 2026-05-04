#include "lancet/base/eta_timer.h"

#include "lancet/base/timer.h"
#include "lancet/base/types.h"

#include "absl/time/clock.h"
#include "absl/time/time.h"

namespace lancet::base {

EtaTimer::EtaTimer(usize const num_iterations) : EtaTimer(num_iterations, &absl::Now) {}

EtaTimer::EtaTimer(usize const num_iterations, Timer::ClockFn const clock)
    : mProgressTimer(clock), mNumTotal(num_iterations) {}

void EtaTimer::Increment() {
  mNumDone++;
  mRunStats.Add(absl::ToInt64Nanoseconds(mProgressTimer.Runtime()));
  mProgressTimer.Reset();
}

auto EtaTimer::EstimatedEta() const -> absl::Duration {
  auto const estimated_ns_remaining = static_cast<f64>(mNumTotal - mNumDone) * mRunStats.Mean();
  return absl::Nanoseconds(estimated_ns_remaining);
}

auto EtaTimer::RatePerSecond() const -> f64 {
  // Dimensional analysis:
  //   NS_PER_SECOND = nanosecs / second
  //   mRunStats.Mean() = nanosecs / window
  //   Result = (nanosecs/second) / (nanosecs/window)
  //         ===> windows/second
  static constexpr f64 NS_PER_SECOND = 1e9;
  return NS_PER_SECOND / mRunStats.Mean();
}

}  // namespace lancet::base
