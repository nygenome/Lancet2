#ifndef SRC_LANCET_BASE_ASSERT_H_
#define SRC_LANCET_BASE_ASSERT_H_

#include "spdlog/fmt/bundled/format.h"
#ifdef LANCET_DEBUG_MODE
#include <source_location>
#include <stdexcept>
#endif

#ifdef LANCET_DEBUG_MODE
inline void ThrowIfAssertFail(
    bool condition, std::source_location const location = std::source_location::current()) {
  if (!condition) {
    throw std::runtime_error(fmt::format("assertion failed: {}:{}:{} - `{}`", location.file_name(),
                                         location.line(), location.column(),
                                         location.function_name()));
  }
}

// Macro form is required so LANCET_ASSERT compiles to a no-op in release without evaluating its
// argument; an inline function would still evaluate the expression.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define LANCET_ASSERT(condition) ThrowIfAssertFail(condition);
#else
// Macro form is required so LANCET_ASSERT compiles to a no-op in release without evaluating its
// argument; an inline function would still evaluate the expression.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define LANCET_ASSERT(e) ((void)0);
#endif

#endif  // SRC_LANCET_BASE_ASSERT_H_
