#pragma once

#include <string>

#include "absl/time/clock.h"
#include "absl/time/time.h"

namespace lancet {
[[nodiscard]] inline auto Humanized(const absl::Duration& d) -> std::string { return absl::FormatDuration(d); }

class Timer {
 public:
  Timer() : startTime(absl::Now()) {}

  [[nodiscard]] auto Runtime() -> absl::Duration {
    if (runComplete) return prevElapsed;

    const auto elapsed = absl::Now() - startTime;
    runComplete = true;
    prevElapsed = elapsed;
    return elapsed;
  }

  [[nodiscard]] auto HumanRuntime() -> std::string { return Humanized(Runtime()); }

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
