#ifndef SRC_LANCET_BASE_CRASH_HANDLER_H_
#define SRC_LANCET_BASE_CRASH_HANDLER_H_

#include "lancet/base/types.h"

namespace lancet::base {

/// Install a process-wide SIGSEGV/SIGABRT handler that prints diagnostic info
/// to stderr on crash: faulting address, signal code, thread ID, instruction
/// pointer (RIP), and a best-effort stack backtrace.
///
/// Call once from main() before spawning worker threads. On platforms without
/// POSIX signals (non-Linux, non-macOS) this is a no-op.
void InstallCrashHandler();

// ============================================================================
// Per-thread crash context — cooperatively populated by worker threads.
//
// AsyncWorker::Process() updates the calling thread's crash context before
// each ProcessWindow() call.  When the crash handler fires, it reads all
// registered thread slots to report what each worker was doing at crash time.
//
// The design respects the layer boundary:
//   lancet::base  (layer 1) — owns the slot storage + signal-safe read path.
//   lancet::core  (layer 5) — the only layer that writes to the slots.
//
// Thread lifecycle:
//   1. Thread calls RegisterThreadSlot() once at startup → gets a slot index.
//   2. Loop: SetSlotWindowInfo() before ProcessWindow(), ClearSlotWindowInfo()
//      after.  If a crash occurs mid-processing the slot holds the active window.
//   3. Thread calls UnregisterThreadSlot() before exit.
//
// Implementation: fixed-size atomic slot array — no heap allocation, no locks,
// fully async-signal-safe to read from the crash handler.
// ============================================================================

/// Maximum concurrent worker threads that can register crash context.
/// Pipeline typically uses 1–128 worker threads.  Excess threads beyond
/// this limit silently skip registration (crash context will say "unregistered").
static constexpr usize MAX_CRASH_SLOTS = 256;

/// Maximum length of the human-readable region string stored per slot.
/// Truncated to this length if longer (e.g., "chr4:12345-67890" = 18 chars).
static constexpr usize MAX_REGION_LEN = 128;

/// Opaque slot index returned by RegisterThreadSlot.  A value of SIZE_MAX
/// indicates registration failure (all slots occupied).
using CrashSlotIdx = usize;
static constexpr CrashSlotIdx INVALID_CRASH_SLOT = SIZE_MAX;

/// Allocate a crash context slot for the calling thread.
/// Call once per worker thread, before the processing loop begins.
/// Returns the slot index (used in all subsequent calls) or INVALID_CRASH_SLOT.
auto RegisterThreadSlot() -> CrashSlotIdx;

/// Release a crash context slot.  Call when the worker thread exits.
void UnregisterThreadSlot(CrashSlotIdx slot_idx);

/// Update a slot with the window currently being processed.
/// genome_idx: Window::GenomeIndex()  (globally unique window ordinal)
/// region_str: Window::ToSamtoolsRegion()  (e.g., "chr4:12345-67890")
void SetSlotWindowInfo(CrashSlotIdx slot_idx, u64 genome_idx, char const* region_str);

/// Clear a slot after ProcessWindow() completes (marks thread as idle).
void ClearSlotWindowInfo(CrashSlotIdx slot_idx);

/// RAII wrapper around `RegisterThreadSlot` / `UnregisterThreadSlot`.
/// Construct on a worker thread to claim a crash-context slot; destruct
/// (typically at thread exit) to release it. The free-function API
/// continues to exist for call sites that already manage the slot lifetime
/// directly; new call sites should prefer this scope guard so a return
/// path or thrown exception can't leak the slot.
class CrashContextScope {
 public:
  CrashContextScope() : mSlotIdx(RegisterThreadSlot()) {}

  ~CrashContextScope() {
    if (mSlotIdx != INVALID_CRASH_SLOT) {
      UnregisterThreadSlot(mSlotIdx);
    }
  }

  CrashContextScope(CrashContextScope const&) = delete;
  auto operator=(CrashContextScope const&) -> CrashContextScope& = delete;
  CrashContextScope(CrashContextScope&&) noexcept = delete;
  auto operator=(CrashContextScope&&) noexcept -> CrashContextScope& = delete;

  /// The slot index assigned to the calling thread, or INVALID_CRASH_SLOT
  /// if all slots were occupied at construction time. Callers pass this
  /// index to SetSlotWindowInfo / ClearSlotWindowInfo.
  [[nodiscard]] auto SlotIdx() const noexcept -> CrashSlotIdx { return mSlotIdx; }

  /// True iff a slot was successfully claimed.
  [[nodiscard]] auto IsValid() const noexcept -> bool { return mSlotIdx != INVALID_CRASH_SLOT; }

 private:
  CrashSlotIdx mSlotIdx;
};

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_CRASH_HANDLER_H_
