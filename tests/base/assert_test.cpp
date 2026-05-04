#include "lancet/base/assert.h"

#include "catch_amalgamated.hpp"

#include <stdexcept>

// CRITICAL: This test file verifies the THROW-on-false debug behaviour of
// `LANCET_ASSERT`. The release no-op behavior — `LANCET_ASSERT` expanding to
// `((void)0)` when `LANCET_DEBUG_MODE` is undefined — is NOT verified here
// because the project's test target unconditionally defines
// `LANCET_DEBUG_MODE` (see `tests/CMakeLists.txt`'s
// `target_compile_definitions(TestLancet2 PRIVATE LANCET_DEBUG_MODE)`).
// Verifying the release no-op would require a separate translation unit
// without the macro defined; that's outside the unit-test target and
// covered by inspection of `assert.h` itself.

namespace lancet::base::tests {

TEST_CASE("LANCET_ASSERT(true) is a no-op under LANCET_DEBUG_MODE",
          "[lancet][base][LancetAssert]") {
  // The non-throwing branch must execute cleanly. `ThrowIfAssertFail(true)`
  // returns without throwing; the macro expands to a single statement.
  CHECK_NOTHROW([]() { LANCET_ASSERT(true) }());
}

TEST_CASE("LANCET_ASSERT(false) throws std::runtime_error under LANCET_DEBUG_MODE",
          "[lancet][base][LancetAssert]") {
  // The failing branch throws std::runtime_error with the source location
  // embedded in the message. We assert the exception type — the exact
  // message format is implementation-detail and would over-couple this test
  // to the formatter.
  CHECK_THROWS_AS([]() { LANCET_ASSERT(false) }(), std::runtime_error);
}

TEST_CASE("LANCET_ASSERT does not evaluate side effects more than once",
          "[lancet][base][LancetAssert]") {
  // The macro expands to `ThrowIfAssertFail(condition)` — one evaluation of
  // the condition expression. Smoke-test that an expression with a side
  // effect runs exactly once, i.e. the macro doesn't re-evaluate via some
  // do-while idiom.
  int call_count = 0;
  auto const tick_and_check = [&]() {
    ++call_count;
    return true;
  };
  LANCET_ASSERT(tick_and_check())
  CHECK(call_count == 1);
}

}  // namespace lancet::base::tests
