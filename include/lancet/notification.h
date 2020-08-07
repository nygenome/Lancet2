#pragma once

#include <stdexcept>
#include <string>

#include "absl/synchronization/notification.h"
#include "absl/time/time.h"
#include "lancet/timer.h"

namespace lancet {
template <typename T>
class Notification : private absl::Notification {
 public:
  Notification() = default;

  [[nodiscard]] auto IsNotDone() const -> bool { return !IsDone(); }
  [[nodiscard]] auto IsDone() const -> bool { return HasBeenNotified(); }
  [[nodiscard]] auto HumanRuntime() const -> std::string { return Humanized(duration); }

  [[nodiscard]] auto WaitForResult() const -> T {
    WaitForNotification();
    return result;
  }

  void SetResult(const T& val, const absl::Duration& runtime) {
    result = val;
    duration = runtime;
    if (!HasBeenNotified()) Notify();
  }

 private:
  T result;
  absl::Duration duration;
};
}  // namespace lancet
