#include "lancet/base/memory.h"

#include "catch_amalgamated.hpp"

#include <string>

namespace lancet::base::tests {

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PeakRssBytes — process-cumulative peak RSS                              ║
// ║                                                                          ║
// ║  CRITICAL: PeakRssBytes is process-cumulative via                        ║
// ║  getrusage(RUSAGE_SELF) and reflects whatever earlier tests in the       ║
// ║  same suite allocated. Any "after a specific allocation" assertion       ║
// ║  would be order-dependent — a regression-safe assertion is bounded to    ║
// ║  "nonzero" (the kernel always reports SOME max RSS for a running         ║
// ║  process).                                                               ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("PeakRssBytes returns a nonzero value for the running process",
          "[lancet][base][PeakRssBytes]") {
  // The running test process has been alive for at least the time it took
  // Catch2 to register every TEST_CASE; getrusage(RUSAGE_SELF) reports the
  // peak resident-set-size and that must be nonzero. The exact value is
  // determined by whatever the suite has allocated up to this point, so we
  // do not assert any specific magnitude.
  auto const bytes = PeakRssBytes();
  CHECK(bytes > 0);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  FormatPeakRss — pure human-readable formatter                           ║
// ║                                                                          ║
// ║  FormatPeakRss reads PeakRssBytes() internally, so the exact unit it     ║
// ║  picks (MB vs GB) is process-state-dependent. We assert structural       ║
// ║  invariants that hold regardless: the formatted string is non-empty      ║
// ║  and ends in a recognized unit suffix.                                   ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("FormatPeakRss produces a non-empty string with a unit suffix",
          "[lancet][base][FormatPeakRss]") {
  auto const formatted = FormatPeakRss();
  CHECK_FALSE(formatted.empty());

  // FormatPeakRss returns "%.2f GB" if peak ≥ 1 GiB, otherwise "%.1f MB".
  // Either suffix must be present at the very end.
  auto const ends_with_gb = formatted.ends_with(" GB");
  auto const ends_with_mb = formatted.ends_with(" MB");
  INFO("formatted=\"" << formatted << "\"");
  CHECK((ends_with_gb || ends_with_mb));
}

}  // namespace lancet::base::tests
