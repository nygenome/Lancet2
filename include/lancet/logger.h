#pragma once

#include <iostream>
#include <utility>

#include "absl/strings/str_format.h"
#include "absl/synchronization/mutex.h"

namespace lancet {
/// Returns RFC3339_sec formatted local time
auto RFC3339Time() -> std::string;

static absl::Mutex logLock;  // NOLINT

/// Prints the message formatted with the provided args and prefixed with RFC3339 timestamp
template <typename... Args>
inline auto InfoLog(const absl::FormatSpec<Args...>& format, Args&&... args) -> void {
  absl::MutexLock guard(&logLock);
  std::cerr << absl::StreamFormat("%s - INFO - %s\n", RFC3339Time(),
                                  absl::StrFormat(format, std::forward<Args>(args)...));
}

template <typename... Args>
inline auto WarnLog(const absl::FormatSpec<Args...>& format, Args&&... args) -> void {
  absl::MutexLock guard(&logLock);
  std::cerr << absl::StreamFormat("%s - WARN - %s\n", RFC3339Time(),
                                  absl::StrFormat(format, std::forward<Args>(args)...));
}

template <typename... Args>
inline auto FatalLog(const absl::FormatSpec<Args...>& format, Args&&... args) -> void {
  absl::MutexLock guard(&logLock);
  std::cerr << absl::StreamFormat("%s - FATA - %s\n", RFC3339Time(),
                                  absl::StrFormat(format, std::forward<Args>(args)...));
}

#if !defined(NDEBUG)
template <typename... Args>
inline auto DebugLog(const absl::FormatSpec<Args...>& format, Args&&... args) -> void {
  absl::MutexLock guard(&logLock);
  std::cerr << absl::StreamFormat("%s - DBUG - %s\n", RFC3339Time(),
                                  absl::StrFormat(format, std::forward<Args>(args)...));
}
#else
#define DebugLog(format, ...) (void)0  // NOLINT
#endif
}  // namespace lancet
