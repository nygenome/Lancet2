#pragma once

#include <string>

#include "absl/time/clock.h"
#include "absl/time/time.h"

namespace lancet {
class Timer {
 public:
  Timer() : startTime(absl::Now()) {}

  auto Runtime() -> absl::Duration {
    if (runComplete) return prevElapsed;

    const auto elapsed = absl::Now() - startTime;
    runComplete = true;
    prevElapsed = elapsed;
    return elapsed;
  }

  auto HumanRuntime() -> std::string { return absl::FormatDuration(Runtime()); }

  void Reset() {
    runComplete = false;
    startTime = absl::Now();
  }

 private:
  absl::Time startTime;
  absl::Duration prevElapsed;
  bool runComplete = false;
};
}  // namespace lancet
