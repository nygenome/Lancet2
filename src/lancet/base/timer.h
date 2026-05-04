#ifndef SRC_LANCET_BASE_TIMER_H_
#define SRC_LANCET_BASE_TIMER_H_

#include "absl/time/clock.h"
#include "absl/time/time.h"

#include <string>

namespace lancet::base {

class Timer {
 public:
  // Function-pointer clock seam. The default ctor wires `&absl::Now`, so
  // production callers see no behavior change. A test ctor injecting a fixture
  // clock makes every "now" read inside Runtime/HumanRuntime/Reset deterministic.
  // A start-time-only seam is insufficient — those methods read the current
  // time during their work, not at construction.
  using ClockFn = absl::Time (*)();

  Timer() : Timer(&absl::Now) {}
  explicit Timer(ClockFn clock) : mClock(clock), mStartTime(clock()) {}

  [[nodiscard]] auto Runtime() -> absl::Duration { return mClock() - mStartTime; }
  [[nodiscard]] auto HumanRuntime() -> std::string { return absl::FormatDuration(Runtime()); }

  void Reset() { mStartTime = mClock(); }

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  ClockFn mClock;
  absl::Time mStartTime;
};

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_TIMER_H_
