#pragma once

#ifndef NDEBUG
#include <stdexcept>

#include "absl/strings/str_format.h"
// NOLINTNEXTLINE
#define LANCET_ASSERT(condition) \
  if (!(condition)) throw std::runtime_error(absl::StrFormat("assertion failed: %s:%d", __FILE__, __LINE__));
#else

#define LANCET_ASSERT(condition) (void)0  // NOLINT
#endif
