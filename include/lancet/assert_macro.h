#pragma once

#ifndef NDEBUG
#include <stdexcept>
// NOLINTNEXTLINE
#define LANCET_ASSERT(condition) \
  if (!(condition)) throw std::runtime_error("assertion failed!");
#else

#define LANCET_ASSERT(condition) (void)0  // NOLINT
#endif
