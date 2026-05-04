#include "lancet/base/timer.h"

#include "absl/time/time.h"
#include "catch_amalgamated.hpp"

#include <atomic>
#include <string>

namespace lancet::base::tests {

namespace {

// Mock clock — every call to Now() returns the value at *gActiveClockTick*.
// Tests advance the tick to script the timer's view of the world without
// touching the real wall clock. The pointer plumbing lets us use a single
// function-pointer signature (matching `Timer::ClockFn`) while letting each
// TEST_CASE own its own ticks. Atomic for theoretical thread-safety; not
// strictly required since each TEST_CASE runs single-threaded.
//
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::atomic<absl::Time const*> gActiveClockTick{nullptr};

auto MockNow() -> absl::Time {
  // Each TEST_CASE installs its own clock-tick pointer. Reading nullptr is
  // a programmer error; surfacing it via a sentinel time keeps the failure
  // diagnosable rather than crashing on a deref.
  auto const* const tick = gActiveClockTick.load(std::memory_order_acquire);
  return tick != nullptr ? *tick : absl::UnixEpoch();
}

}  // namespace

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  Timer — clock-seam tests                                                ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("Timer Runtime returns the duration since construction under the injected clock",
          "[lancet][base][Timer]") {
  // Start the clock at a fixed instant; advance by exactly 5s before the
  // first Runtime() read. The clock seam guarantees Runtime() observes the
  // advanced tick — without the seam, this test would race the wall clock.
  auto current = absl::FromUnixSeconds(1'000'000);
  gActiveClockTick.store(&current, std::memory_order_release);
  Timer timer(&MockNow);

  current = absl::FromUnixSeconds(1'000'005);
  CHECK(timer.Runtime() == absl::Seconds(5));

  current = absl::FromUnixSeconds(1'000'042);
  CHECK(timer.Runtime() == absl::Seconds(42));

  gActiveClockTick.store(nullptr, std::memory_order_release);
}

TEST_CASE("Timer Reset rebases the start to the current clock tick", "[lancet][base][Timer]") {
  auto current = absl::FromUnixSeconds(1'000'000);
  gActiveClockTick.store(&current, std::memory_order_release);
  Timer timer(&MockNow);

  current = absl::FromUnixSeconds(1'000'010);
  CHECK(timer.Runtime() == absl::Seconds(10));

  // Reset rebases mStartTime to the current tick. Subsequent Runtime() is
  // measured from this new origin, not the original construction point.
  timer.Reset();
  CHECK(timer.Runtime() == absl::ZeroDuration());

  current = absl::FromUnixSeconds(1'000'013);
  CHECK(timer.Runtime() == absl::Seconds(3));

  gActiveClockTick.store(nullptr, std::memory_order_release);
}

TEST_CASE("Timer HumanRuntime renders the duration via absl::FormatDuration",
          "[lancet][base][Timer]") {
  // HumanRuntime delegates to absl::FormatDuration, so we verify only that
  // the produced string is the formatted form of the controllable duration —
  // not the exact format string (which is absl::FormatDuration's contract,
  // not Timer's).
  auto current = absl::FromUnixSeconds(1'000'000);
  gActiveClockTick.store(&current, std::memory_order_release);
  Timer timer(&MockNow);

  current = absl::FromUnixSeconds(1'000'002);  // +2s
  auto const text = timer.HumanRuntime();
  CHECK(text == absl::FormatDuration(absl::Seconds(2)));

  gActiveClockTick.store(nullptr, std::memory_order_release);
}

TEST_CASE("Timer default ctor wires absl::Now and produces non-negative runtime",
          "[lancet][base][Timer]") {
  // Production callers see no behavior change. The default ctor's "now" is
  // not deterministic, so we only assert the post-condition Runtime() >= 0.
  Timer timer;
  CHECK(timer.Runtime() >= absl::ZeroDuration());
  CHECK_FALSE(timer.HumanRuntime().empty());
}

}  // namespace lancet::base::tests
