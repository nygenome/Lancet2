#ifndef SRC_LANCET_BASE_MEMORY_H_
#define SRC_LANCET_BASE_MEMORY_H_

#include "lancet/base/types.h"

#include "spdlog/fmt/bundled/format.h"

// POSIX header — quote-included and `extern "C"`-wrapped per the project
// convention for non-stdlib C headers (cpp_style.md, mirrors lancet/hts/).
// IWYU pragma: no_include <sys/resource.h>
extern "C" {
#include "sys/resource.h"  // IWYU pragma: keep
}

#include <string>

namespace lancet::base {

// ============================================================================
// PeakRssBytes — peak resident set size (max RSS) of the current process.
//
// Uses getrusage(RUSAGE_SELF) which reports ru_maxrss in:
//   Linux:  kilobytes (× 1024 to get bytes)
//   macOS:  bytes     (no conversion needed)
//
// Returns 0 on failure (getrusage returns -1).
// ============================================================================
[[nodiscard]] inline auto PeakRssBytes() -> usize {
  struct rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;

  // POSIX `struct rusage` declares `ru_maxrss` inside a tagged union on some libc implementations
  // for cross-platform layout compatibility; direct member access is the documented contract.
  // NOLINTBEGIN(cppcoreguidelines-pro-type-union-access)
#ifdef __APPLE__
  return static_cast<usize>(usage.ru_maxrss);  // macOS: already in bytes
#else
  return static_cast<usize>(usage.ru_maxrss) * 1024;  // Linux: kilobytes → bytes
#endif
  // NOLINTEND(cppcoreguidelines-pro-type-union-access)
}

// ============================================================================
// FormatPeakRss — human-readable peak RSS string (e.g., "1.23 GB").
// ============================================================================
[[nodiscard]] inline auto FormatPeakRss() -> std::string {
  auto const bytes = PeakRssBytes();
  static constexpr f64 BYTES_PER_GB = 1024.0 * 1024.0 * 1024.0;
  static constexpr f64 BYTES_PER_MB = 1024.0 * 1024.0;

  if (auto const mem_in_gb = static_cast<f64>(bytes) / BYTES_PER_GB; mem_in_gb >= 1.0) {
    return fmt::format("{:.2f} GB", mem_in_gb);
  }

  auto const mem_in_mb = static_cast<f64>(bytes) / BYTES_PER_MB;
  return fmt::format("{:.1f} MB", mem_in_mb);
}

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_MEMORY_H_
