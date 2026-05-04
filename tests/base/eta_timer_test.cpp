#include "lancet/base/eta_timer.h"

#include "lancet/base/types.h"

#include "absl/time/time.h"
#include "catch_amalgamated.hpp"

#include <atomic>

namespace lancet::base::tests {

namespace {

// Same mock-clock pattern as timer_test.cpp: a globally-pointed tick that
// the function-pointer clock dereferences on every call. EtaTimer threads
// the clock through its internal Timer, so this mock makes both Increment()
// and EstimatedEta() deterministic.
//
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::atomic<absl::Time const*> gActiveClockTick{nullptr};

auto MockNow() -> absl::Time {
  auto const* const tick = gActiveClockTick.load(std::memory_order_acquire);
  return tick != nullptr ? *tick : absl::UnixEpoch();
}

}  // namespace

TEST_CASE("EtaTimer Increment advances the per-iteration runtime under the injected clock",
          "[lancet][base][EtaTimer]") {
  // 10 iterations total. Each iteration consumes 1s on the mock clock. After
  // 3 iterations, RatePerSecond should approach 1.0 and the ETA should match
  // hand-computed (7 iterations * 1s).
  auto current = absl::FromUnixSeconds(1'000'000);
  gActiveClockTick.store(&current, std::memory_order_release);
  EtaTimer eta(/*num_iterations=*/10, &MockNow);

  for (i32 step = 1; step <= 3; ++step) {
    current = absl::FromUnixSeconds(1'000'000 + step);
    eta.Increment();
  }

  CHECK(eta.RatePerSecond() == Catch::Approx(1.0).margin(1e-9));
  CHECK(eta.EstimatedEta() == absl::Seconds(7));

  gActiveClockTick.store(nullptr, std::memory_order_release);
}

TEST_CASE("EtaTimer EstimatedEta is monotone-decreasing as work completes",
          "[lancet][base][EtaTimer]") {
  // ETA after k iterations is (n-k)*mean_step. Mean stays the same when
  // every step takes the same time, so ETA decreases linearly toward 0.
  // This guards against a regression that resets the running mean on each
  // Increment (which would freeze the ETA) or that confuses (n-k) with k.
  auto current = absl::FromUnixSeconds(1'000'000);
  gActiveClockTick.store(&current, std::memory_order_release);
  EtaTimer eta(/*num_iterations=*/5, &MockNow);

  current = absl::FromUnixSeconds(1'000'001);
  eta.Increment();
  auto const eta_after_1 = eta.EstimatedEta();

  current = absl::FromUnixSeconds(1'000'002);
  eta.Increment();
  auto const eta_after_2 = eta.EstimatedEta();

  current = absl::FromUnixSeconds(1'000'003);
  eta.Increment();
  auto const eta_after_3 = eta.EstimatedEta();

  CHECK(eta_after_1 > eta_after_2);
  CHECK(eta_after_2 > eta_after_3);
  CHECK(eta_after_3 == absl::Seconds(2));  // 5 - 3 remaining * 1s/step

  gActiveClockTick.store(nullptr, std::memory_order_release);
}

TEST_CASE("EtaTimer RatePerSecond matches hand-computed reciprocal of mean step time",
          "[lancet][base][EtaTimer]") {
  // Two 2-second steps → mean step = 2s → rate = 0.5 windows/s. The
  // dimensional analysis in the EtaTimer source documents this calculation;
  // the test pins it to a known reference.
  auto current = absl::FromUnixSeconds(1'000'000);
  gActiveClockTick.store(&current, std::memory_order_release);
  EtaTimer eta(/*num_iterations=*/4, &MockNow);

  current = absl::FromUnixSeconds(1'000'002);
  eta.Increment();
  current = absl::FromUnixSeconds(1'000'004);
  eta.Increment();

  CHECK(eta.RatePerSecond() == Catch::Approx(0.5).margin(1e-9));

  gActiveClockTick.store(nullptr, std::memory_order_release);
}

TEST_CASE("EtaTimer default ctor wires absl::Now without affecting Increment correctness",
          "[lancet][base][EtaTimer]") {
  // Production path uses the default ctor. The wall clock makes the rate
  // non-deterministic, so we only assert that Increment progresses state
  // and EstimatedEta becomes a finite (non-infinite) absl::Duration after
  // at least one increment.
  EtaTimer eta(/*num_iterations=*/3);
  eta.Increment();
  auto const computed_eta = eta.EstimatedEta();
  CHECK(computed_eta >= absl::ZeroDuration());
  CHECK(computed_eta != absl::InfiniteDuration());
}

}  // namespace lancet::base::tests
