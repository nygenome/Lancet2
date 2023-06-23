#ifndef SRC_LANCET_BASE_TIMER_H_
#define SRC_LANCET_BASE_TIMER_H_

#include <string>

#include "absl/time/clock.h"
#include "absl/time/time.h"

class Timer {
 public:
  Timer() : mStartTime(absl::Now()) {}

  [[nodiscard]] auto Runtime() -> absl::Duration { return absl::Now() - mStartTime; }
  [[nodiscard]] auto HumanRuntime() -> std::string { return absl::FormatDuration(Runtime()); }

  void Reset() { mStartTime = absl::Now(); }

 private:
  absl::Time mStartTime;
};

#endif  // SRC_LANCET_BASE_TIMER_H_
