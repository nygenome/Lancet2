#include "lancet/base/eta_timer.h"

#include "lancet/base/types.h"

#include "absl/time/time.h"

#include <cmath>

namespace lancet::base {

EtaTimer::EtaTimer(usize const num_iterations) : mNumTotal(num_iterations) {}

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
  static constexpr f64 NS_TO_SECS = 1e-9;
  static constexpr f64 WINDOWS_PER_SECOND_CONVERTER = -1.0;
  return std::pow(mRunStats.Mean() * NS_TO_SECS, WINDOWS_PER_SECOND_CONVERTER);
}

}  // namespace lancet::base
