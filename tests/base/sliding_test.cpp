#include "lancet/base/sliding.h"

#include "catch_amalgamated.hpp"

#include <string>
#include <string_view>

namespace lancet::base::tests {

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  SlidingView — k-mer enumeration                                         ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("SlidingView yields seq.length() − window + 1 windows when window < seq",
          "[lancet][base][SlidingView]") {
  // Canonical case. For a sequence of length n and window size k, the number
  // of valid windows is (n - k + 1). Each window is an exclusive slice of
  // length k. The implementation guarantees this via an internal LANCET_ASSERT.
  auto const windows = SlidingView("ACGTACGT", 3);
  REQUIRE(windows.size() == 6);  // 8 − 3 + 1
  CHECK(windows[0] == "ACG");
  CHECK(windows[1] == "CGT");
  CHECK(windows[2] == "GTA");
  CHECK(windows[3] == "TAC");
  CHECK(windows[4] == "ACG");
  CHECK(windows[5] == "CGT");
}

TEST_CASE("SlidingView yields exactly one window when window == seq",
          "[lancet][base][SlidingView]") {
  // Smallest non-empty result. A sequence of length k slid with window k
  // yields one window equal to the whole sequence. Off-by-one regressions
  // (returning 0 or 2 windows) trip here.
  auto const windows = SlidingView("ACGT", 4);
  REQUIRE(windows.size() == 1);
  CHECK(windows[0] == "ACGT");
}

TEST_CASE("SlidingView yields an empty result when window > seq", "[lancet][base][SlidingView]") {
  // A window larger than the sequence cannot fit; the implementation returns
  // a zero-size FixedArray rather than throwing or under-flowing the size_t
  // arithmetic (`seq.length() - window` would otherwise wrap to a huge value).
  auto const windows = SlidingView("ACG", 5);
  CHECK(windows.empty());
}

TEST_CASE("SlidingView yields an empty result when seq is empty", "[lancet][base][SlidingView]") {
  // Empty input. The window-greater-than-seq guard catches this case (any
  // window > 0 exceeds an empty seq); a window of 0 would be ill-formed
  // (zero-length slice) and isn't a documented input.
  auto const windows = SlidingView("", 3);
  CHECK(windows.empty());
}

TEST_CASE("SlidingView handles window=1 (single-base windows)", "[lancet][base][SlidingView]") {
  // Edge case: window-of-one yields one window per character. Exercises the
  // path where the offset range is the full sequence length.
  auto const windows = SlidingView("ACGT", 1);
  REQUIRE(windows.size() == 4);
  CHECK(windows[0] == "A");
  CHECK(windows[1] == "C");
  CHECK(windows[2] == "G");
  CHECK(windows[3] == "T");
}

}  // namespace lancet::base::tests
