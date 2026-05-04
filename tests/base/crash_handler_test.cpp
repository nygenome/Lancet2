#include "lancet/base/crash_handler.h"

#include "catch_amalgamated.hpp"

#include <string>

namespace lancet::base::tests {

// CRITICAL: The signal-handler portion of crash_handler.{h,cpp} —
// `InstallCrashHandler()` and the `CrashHandler` body that runs on SIGSEGV /
// SIGABRT — is integration-only and not exercised by these unit tests.
// Installing process-wide signal handlers, formatting on an alternate signal
// stack, and re-raising for a core dump are all global-state behaviors that
// are tested end-to-end via the project's e2e crash harness, not here. These
// tests cover only the slot-management primitives and the new
// `CrashContextScope` RAII guard.

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  CrashContextScope — RAII guard                                          ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("CrashContextScope claims a valid slot index by default",
          "[lancet][base][CrashContextScope]") {
  // Fresh scope on a previously-untouched slot table claims a valid index.
  // After destruction the slot is released; subsequent scopes can reclaim
  // any freed slot, exercising the register/unregister round-trip.
  CrashSlotIdx claimed_idx = INVALID_CRASH_SLOT;
  {
    CrashContextScope const scope;
    CHECK(scope.IsValid());
    CHECK(scope.SlotIdx() != INVALID_CRASH_SLOT);
    CHECK(scope.SlotIdx() < MAX_CRASH_SLOTS);
    claimed_idx = scope.SlotIdx();
  }
  // After destruction, claiming a new scope reuses the same slot (slots
  // are claimed first-fit from the front of the array).
  CrashContextScope const reclaimed;
  CHECK(reclaimed.IsValid());
  CHECK(reclaimed.SlotIdx() == claimed_idx);
}

TEST_CASE("CrashContextScope nesting allocates distinct slots",
          "[lancet][base][CrashContextScope]") {
  // Two simultaneously-live scopes occupy two distinct slots. This is the
  // multi-thread case in miniature — each "thread" claims its own slot.
  CrashContextScope const outer;
  CrashContextScope const inner;
  REQUIRE(outer.IsValid());
  REQUIRE(inner.IsValid());
  CHECK(outer.SlotIdx() != inner.SlotIdx());
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  Slot-info round-trip through Set/ClearSlotWindowInfo                    ║
// ║                                                                          ║
// ║  The slot's mGenomeIdx and mRegion fields are populated cooperatively by ║
// ║  the worker thread and read async-signal-safely by the crash handler.    ║
// ║  These tests validate the set/clear primitives without invoking the      ║
// ║  signal handler — the handler reads the same fields through              ║
// ║  DumpAllThreadContexts, which is integration-only.                       ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("SetSlotWindowInfo / ClearSlotWindowInfo are no-ops on INVALID_CRASH_SLOT",
          "[lancet][base][SetSlotWindowInfo]") {
  // Defensive: passing INVALID_CRASH_SLOT must not crash. This protects
  // workers that exhaust the slot table — they keep running with a "no
  // slot available" state, and Set/Clear silently no-op rather than fault.
  CHECK_NOTHROW(SetSlotWindowInfo(INVALID_CRASH_SLOT, 0, "ignored"));
  CHECK_NOTHROW(ClearSlotWindowInfo(INVALID_CRASH_SLOT));
}

TEST_CASE("SetSlotWindowInfo / ClearSlotWindowInfo run without throwing on a live slot",
          "[lancet][base][SetSlotWindowInfo]") {
  // Black-box test: the slot internals are not publicly readable (the
  // signal handler reads them via internal-linkage globals), so we verify
  // only that the public surface accepts the calls and does not throw.
  // The signal-handler path that consumes the data is integration-only.
  CrashContextScope const scope;
  REQUIRE(scope.IsValid());

  CHECK_NOTHROW(SetSlotWindowInfo(scope.SlotIdx(), 42, "chr4:12345-67890"));
  CHECK_NOTHROW(ClearSlotWindowInfo(scope.SlotIdx()));
}

TEST_CASE("SetSlotWindowInfo handles an oversize region string by truncating",
          "[lancet][base][SetSlotWindowInfo]") {
  // The slot stores up to MAX_REGION_LEN bytes including the null terminator.
  // A longer region string must be truncated, not overflow the buffer (the
  // implementation copies up to MAX_REGION_LEN-1 chars and writes the
  // null). Test only that the call doesn't throw — the truncation is
  // observed through the signal handler, which is integration-only.
  CrashContextScope const scope;
  REQUIRE(scope.IsValid());

  // Build a region string longer than MAX_REGION_LEN (128 bytes).
  std::string oversize_region;
  oversize_region.reserve(MAX_REGION_LEN + 32);
  while (oversize_region.size() < MAX_REGION_LEN + 32) {
    oversize_region.append("0123456789");
  }
  CHECK_NOTHROW(SetSlotWindowInfo(scope.SlotIdx(), 0, oversize_region.c_str()));
  CHECK_NOTHROW(ClearSlotWindowInfo(scope.SlotIdx()));
}

TEST_CASE("SetSlotWindowInfo handles a null region pointer", "[lancet][base][SetSlotWindowInfo]") {
  // Defensive: nullptr region_str must not deref. The implementation
  // checks for null before iterating the C string. Without this guard a
  // worker that fails to populate the region would crash the process at
  // exactly the moment crash diagnostics are most needed.
  CrashContextScope const scope;
  REQUIRE(scope.IsValid());
  CHECK_NOTHROW(SetSlotWindowInfo(scope.SlotIdx(), 0, nullptr));
  CHECK_NOTHROW(ClearSlotWindowInfo(scope.SlotIdx()));
}

}  // namespace lancet::base::tests
