#pragma once

#include <stdexcept>
#include <string>

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
    if (!errorMsg.empty()) throw std::runtime_error(errorMsg);
    return result;
  }

  void SetResult(const T& val) {
    result = val;
    if (!HasBeenNotified()) Notify();
  }

  void SetErrorMsg(const std::string& msg) {
    errorMsg = msg;
    if (!HasBeenNotified()) Notify();
  }

 private:
  T result;
  std::string errorMsg;
};
}  // namespace lancet
