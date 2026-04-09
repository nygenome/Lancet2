#ifndef SRC_LANCET_BASE_ASSERT_H_
#define SRC_LANCET_BASE_ASSERT_H_

#ifdef LANCET_DEVELOP_MODE
#include "spdlog/fmt/bundled/core.h"

#include <source_location>
#include <stdexcept>
#endif

#ifdef LANCET_DEVELOP_MODE
inline void ThrowIfAssertFail(
    bool condition, std::source_location const location = std::source_location::current()) {
  if (!condition) {
    throw std::runtime_error(fmt::format("assertion failed: {}:{}:{} - `{}`", location.file_name(),
                                         location.line(), location.column(),
                                         location.function_name()));
  }
}

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define LANCET_ASSERT(condition) ThrowIfAssertFail(condition);
#else
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define LANCET_ASSERT(e) ((void)0);
#endif

#endif  // SRC_LANCET_BASE_ASSERT_H_
