#pragma once

#include "absl/synchronization/notification.h"

namespace lancet {
template <typename T>
class Notification : private absl::Notification {
 public:
  Notification() = default;

  [[nodiscard]] auto IsNotDone() const -> bool { return !IsDone(); }
  [[nodiscard]] auto IsDone() const -> bool { return HasBeenNotified(); }

  [[nodiscard]] auto WaitForResult() const -> T {
    WaitForNotification();
    return result;
  }

  void SetResult(const T& val) {
    result = val;
    if (!HasBeenNotified()) Notify();
  }

 private:
  T result;
};
}  // namespace lancet
