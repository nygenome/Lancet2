#include "lancet/base/crash_handler.h"

#include "lancet/base/types.h"

// ============================================================================
// crash_handler.cpp — POSIX signal-based crash diagnostics.
//
// PURPOSE:
//   When the pipeline encounters a fatal signal (SIGSEGV, SIGABRT) that
//   bypasses normal C++ exception handling, this module prints a diagnostic
//   report to stderr showing which worker thread crashed and which genomic
//   window it was processing.  This turns a blank "Segmentation fault" into
//   an actionable single-threaded reproduction command.
//
// HOW IT WORKS (high level):
//
//   ┌──────────────┐   Worker threads register a "crash slot" and update
//   │ Worker N     │   it with the current window before each ProcessWindow()
//   │  ┌─────────┐ │   call.  If a crash occurs mid-processing, the signal
//   │  │ Slot[N] │ │   handler reads all slots and prints their state.
//   │  └─────────┘ │
//   └──────┬───────┘
//          │ SIGSEGV / SIGABRT
//          ▼
//   ┌──────────────┐   The signal handler runs on a dedicated 64KB alternate
//   │ CrashHandler │   stack (so it works even on stack overflow), prints
//   │  - signal    │   diagnostic info, then re-raises the signal to produce
//   │  - backtrace │   a core dump via the default handler.
//   │  - contexts  │
//   └──────────────┘
//
// IMPORTANT CONSTRAINT — async-signal safety:
//   Code inside a signal handler must NOT call malloc, printf, or any function
//   that takes a lock.  This is because the crash may have occurred while the
//   thread was holding a lock or inside malloc — calling those functions again
//   would deadlock.  All output in this file uses write() (which is safe) and
//   manual formatting into stack-allocated buffers (no heap).
//
// FILE ORGANIZATION:
//   §1  Signal-safe output helpers    — WriteStr, WriteHex, WriteInt, WriteU64
//   §2  Per-thread crash context      — CrashSlot array + DumpAllThreadContexts
//   §3  All-thread instruction pointers — DumpAllThreadIPs (Linux only, via /proc)
//   §4  Signal handler                — CrashHandler (the actual SIGSEGV handler)
//   §5  Public API                    — InstallCrashHandler, Register/Set/Clear
// ============================================================================

#if defined(__linux__) || defined(__APPLE__)

// macOS marks the ucontext routines as deprecated and gates the declarations
// behind _XOPEN_SOURCE.  Define the feature-test macro before any system
// header so that <ucontext.h> (pulled in transitively by <signal.h> as well
// as directly below) exposes the full ucontext_t definition we need to read
// the faulting instruction pointer from the machine context.
#ifdef __APPLE__
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#endif

#include <execinfo.h>
#include <fcntl.h>
#include <pthread.h>
#include <ucontext.h>
#include <unistd.h>

#ifdef __linux__
#include <sys/syscall.h>
#endif

#include <array>
#include <vector>

#include <csignal>
#include <cstring>

namespace {

constexpr int MAX_BACKTRACE_FRAMES = 64;

// Minimum alternate stack size.  We want at least 64KB so the handler has
// room for backtrace() and our formatting buffers.  SIGSTKSZ is checked at
// runtime (glibc >= 2.34 made it a sysconf() call, so it is no longer a
// compile-time constant).
constexpr usize MIN_ALT_STACK_SIZE = 65'536;

auto GetAltStackSize() -> usize {
  auto const sys_size = static_cast<usize>(SIGSTKSZ);
  return sys_size > MIN_ALT_STACK_SIZE ? sys_size : MIN_ALT_STACK_SIZE;
}

// ============================================================================
// §1  Signal-safe output helpers
//
// These functions format and write data using only write() and stack-allocated
// buffers.  They replace printf/fprintf which are NOT safe to call from a
// signal handler (printf calls malloc internally).
//
// write() is one of the few functions guaranteed safe inside a signal handler
// by POSIX (see `man 7 signal-safety`).  clang-tidy doesn't know this, so
// we suppress its warnings on every write() call.
// ============================================================================

void WriteStr(int file_desc, char const* str) {
  usize len = 0;
  while (str[len] != '\0') {
    len++;
  }
  // write() is async-signal-safe per POSIX
  // NOLINTNEXTLINE(bugprone-signal-handler,cert-sig30-c)
  static_cast<void>(write(file_desc, str, len));
}

// Writes "0x" followed by exactly 16 hex digits (zero-padded).
void WriteHex(int file_desc, u64 val) {
  // stack buffer, no heap
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
  char buf[18] = {'0', 'x'};
  for (int idx = 17; idx >= 2; --idx) {
    // index is masked to 0..15 by `val & 0xFU`; bounded access into the 16-char hex digit table.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
    buf[idx] = "0123456789abcdef"[val & 0xFU];
    val >>= 4U;
  }
  // write() is async-signal-safe per POSIX
  // NOLINTNEXTLINE(bugprone-signal-handler,cert-sig30-c)
  static_cast<void>(write(file_desc, buf, sizeof(buf)));
}

// Writes a signed integer in decimal.
void WriteInt(int file_desc, int val) {
  // stack buffer, no heap
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
  char buf[16];
  int pos = 0;
  if (val < 0) {
    buf[pos++] = '-';
    val = -val;
  }

  int const digit_start = pos;
  // do-while is the natural shape of the itoa idiom: at least one digit must be emitted even when
  // val==0.
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-do-while)
  do {
    buf[pos++] = static_cast<char>('0' + (val % 10));
    val /= 10;
  } while (val > 0);

  for (int lo = digit_start, hi = pos - 1; lo < hi; ++lo, --hi) {
    char const tmp = buf[lo];
    buf[lo] = buf[hi];
    buf[hi] = tmp;
  }

  // write() is async-signal-safe per POSIX
  // NOLINTNEXTLINE(bugprone-signal-handler,cert-sig30-c)
  static_cast<void>(write(file_desc, buf, pos));
}

// Writes an unsigned 64-bit integer in decimal.
void WriteU64(int file_desc, u64 val) {
  // stack buffer, no heap
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
  char buf[20];
  int pos = 0;

  // do-while is the natural shape of the itoa idiom: at least one digit must be emitted even when
  // val==0.
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-do-while)
  do {
    buf[pos++] = static_cast<char>('0' + (val % 10));
    val /= 10;
  } while (val > 0);

  for (int lo = 0, hi = pos - 1; lo < hi; ++lo, --hi) {
    char const tmp = buf[lo];
    buf[lo] = buf[hi];
    buf[hi] = tmp;
  }

  // write() is async-signal-safe per POSIX
  // NOLINTNEXTLINE(bugprone-signal-handler,cert-sig30-c)
  static_cast<void>(write(file_desc, buf, pos));
}

// ============================================================================
// §2  Per-thread crash context (the "crash slot" system)
//
// Each worker thread "registers" a slot in a fixed-size global array.  Before
// processing a genomic window, the thread writes the window's index and region
// string into its slot.  If the thread crashes inside ProcessWindow(), the
// signal handler reads all slots and prints which windows were being processed.
//
// Why not just use try/catch?
//   SIGSEGV and SIGABRT are POSIX signals, not C++ exceptions.  A null pointer
//   dereference or heap corruption does NOT trigger try/catch — it triggers the
//   kernel's signal delivery.  The crash slots capture context for these cases.
//   (C++ exceptions like std::out_of_range ARE caught by the try/catch in
//   AsyncWorker::Process — the crash slots are for everything else.)
//
// Thread safety without locks:
//   Each thread writes ONLY to its own slot, claimed via atomic compare-and-swap
//   (CAS — atomically checks "is this slot free?" and claims it in one step).
//   The signal handler only READS slots.  Atomic stores with release ordering
//   (guarantees that all prior writes are visible before the flag is set)
//   ensure the handler sees consistent data.  No mutex is needed.
//
// Memory layout:
//   ┌──────────────────────────────────────────────────────────────────┐
//   │ CrashSlot (192 bytes, aligned to 64-byte CPU cache line so that  │
//   │ different threads' slots don't share the same cache line, which  │
//   │ would cause unnecessary cross-core synchronization overhead)     │
//   ├──────────┬──────────┬────────────┬────────────┬──────────────────┤
//   │ mInUse   │ mActive  │ mGenomeIdx │ mThreadId  │ mRegion[128]     │
//   │ 4B       │ 4B       │ 8B         │ 8B         │ 128B             │
//   │ 0=free   │ 0=idle   │ window     │ pthread_t  │ chr4:12345-67890 │
//   │ 1=claimed│ 1=active │ ordinal    │            │ NUL-terminated   │
//   └──────────┴──────────┴────────────┴────────────┴──────────────────┘
// ============================================================================

struct alignas(64) CrashSlot {
  u64 volatile mGenomeIdx;  // 8B — Window::GenomeIndex() being processed
  pthread_t mThreadId;      // 8B — owning thread ID for crash correlation
  int volatile mInUse;      // 4B — 0 = free, 1 = claimed by a thread
  int volatile mActive;     // 4B — 0 = idle, 1 = processing a window
  // fixed buffer, async-signal-safe
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
  char mRegion[lancet::base::MAX_REGION_LEN];  // 128B — e.g. "chr4:12345-67890"
};

// mutable by design
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::array<CrashSlot, lancet::base::MAX_CRASH_SLOTS> g_crash_slots = {};

// ============================================================================
// DumpAllThreadContexts: print each worker's state at crash time.
// Called from the signal handler — must be async-signal-safe.
// ============================================================================
void DumpAllThreadContexts(int file_desc) {
  WriteStr(file_desc, "\n── Worker Thread Contexts ──────────────────────\n");

  bool found_any = false;
  for (auto const& slot : g_crash_slots) {
    if (slot.mInUse == 0) continue;
    found_any = true;

    WriteStr(file_desc, "  Thread ");
    // pthread_t is a pointer on macOS
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    WriteHex(file_desc, reinterpret_cast<u64>(slot.mThreadId));

    if (slot.mActive != 0) {
      WriteStr(file_desc, " [ACTIVE] window_idx=");
      WriteU64(file_desc, slot.mGenomeIdx);
      WriteStr(file_desc, " region=");
      WriteStr(file_desc, slot.mRegion);
    } else {
      WriteStr(file_desc, " [idle]");
    }
    WriteStr(file_desc, "\n");
  }

  if (!found_any) {
    WriteStr(file_desc, "  (no worker threads registered)\n");
  }
}

// ============================================================================
// §3  All-thread instruction pointer enumeration (Linux only)
//
// On Linux, /proc/self/task/ lists every OS thread in the process.  For each
// thread, /proc/self/task/<tid>/syscall contains the instruction pointer
// (the address of the machine instruction the thread was executing).
// This lets us see what EVERY thread was doing at crash time — not just the
// workers that registered crash slots, but also the main thread, I/O threads,
// and any profiler threads.
//
// The kernel exposes each thread's state in a text file at the path:
//   /proc/self/task/<tid>/syscall
//
// Example contents when the thread is blocked in a system call:
//   "0 0x7f6... 0x1000 0x0 0x0 0x0 0x0 0x7ffd... 0x7f61..."
//    ^                                            ^SP       ^IP
//    syscall number                               stack ptr instruction ptr
//
// The last two hex values are always the stack pointer and instruction
// pointer (the address of the machine instruction the thread was executing).
// If the thread is running in user code (not inside a kernel call), the
// file contains just the word "running".
//
// This section uses the raw getdents64 Linux syscall to list directory
// entries instead of the standard opendir/readdir functions, because
// opendir calls malloc, which is not safe in a signal handler.
// ============================================================================
#ifdef __linux__

// Kernel dirent64 layout — must match the kernel's struct exactly.
// mirrors kernel struct name
// NOLINTNEXTLINE(readability-identifier-naming)
struct linux_dirent64 {
  u64 mIno;                // 8B — inode number
  u64 mOff;                // 8B — offset to next entry
  unsigned short mReclen;  // 2B — total record length
  unsigned char mType;     // 1B — file type (DT_DIR, DT_REG, etc.)
  // variable-length kernel ABI
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
  char mName[1];  // NUL-terminated filename
};

// ============================================================================
// DumpOneThreadIP: read and print the syscall state for a single thread.
// ============================================================================
void DumpOneThreadIP(int file_desc, char const* tid_str) {
  // Build the path: /proc/self/task/<tid>/syscall
  // stack buffer for path, no heap
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
  char path[64] = "/proc/self/task/";
  constexpr usize PREFIX_LEN = 16;
  usize plen = PREFIX_LEN;

  for (usize idx = 0; tid_str[idx] != '\0' && plen < 48; ++idx) {
    path[plen++] = tid_str[idx];
  }
  memcpy(path + plen, "/syscall", 9);

  // stack buffer for procfs read, no heap
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
  char syscall_buf[256] = {};
  int const sfd = open(path, O_RDONLY);
  if (sfd < 0) return;
  auto const nread = read(sfd, syscall_buf, sizeof(syscall_buf) - 1);
  close(sfd);
  if (nread <= 0) return;

  // Truncate at first newline.
  for (long idx = 0; idx < nread; ++idx) {
    if (syscall_buf[idx] == '\n') {
      syscall_buf[idx] = '\0';
      break;
    }
  }

  WriteStr(file_desc, "  tid=");
  WriteStr(file_desc, tid_str);
  WriteStr(file_desc, " syscall: ");
  WriteStr(file_desc, syscall_buf);
  WriteStr(file_desc, "\n");
}

// ============================================================================
// DumpAllThreadIPs: enumerate all OS threads and print their instruction pointers.
// Uses the raw getdents64 syscall to avoid malloc (opendir is not signal-safe).
// ============================================================================
void DumpAllThreadIPs(int file_desc) {
  WriteStr(file_desc, "\n── All Thread IPs (via /proc/self/task) ────────\n");

  // open() is a POSIX variadic function
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg)
  int const dir_fd = open("/proc/self/task", O_RDONLY | O_DIRECTORY);
  if (dir_fd < 0) {
    WriteStr(file_desc, "  (could not open /proc/self/task)\n");
    return;
  }

  // raw buffer for getdents64 syscall
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
  char dirents_buf[4096];

  while (true) {
    // syscall() is a POSIX variadic function
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg)
    auto const nbytes = syscall(SYS_getdents64, dir_fd, dirents_buf, sizeof(dirents_buf));
    if (nbytes <= 0) break;

    for (long offset = 0; offset < nbytes;) {
      // getdents64 returns variable-length records back-to-back; navigation is by byte offset
      // through `entry->mReclen` per the kernel ABI. The reinterpret_cast is the documented form.
      // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-bounds-pointer-arithmetic)
      auto const* entry = reinterpret_cast<linux_dirent64 const*>(dirents_buf + offset);
      offset += entry->mReclen;

      if (entry->mName[0] == '.') continue;
      DumpOneThreadIP(file_desc, entry->mName);
    }
  }

  close(dir_fd);
}
#endif  // __linux__

// ============================================================================
// §4  Signal handler
//
// This is the function that runs when SIGSEGV or SIGABRT is delivered.
// It prints a diagnostic report, then re-raises the signal with the default
// handler to produce a core dump (a snapshot of the process state that can be
// examined later with a debugger).
//
// The handler uses SA_RESETHAND — it auto-resets to the OS default handler
// (SIG_DFL) before entering.  This means a crash INSIDE the handler itself
// will produce a normal core dump instead of infinite recursion.
// ============================================================================

// signature prescribed by POSIX sigaction
// NOLINTNEXTLINE(cert-err33-c,cert-dcl03-c)
void CrashHandler(int sig, siginfo_t* info, void* ctx) {
  int const file_desc = STDERR_FILENO;
  WriteStr(file_desc, "\n============================================\n");
  WriteStr(file_desc, sig == SIGSEGV ? "*** SIGSEGV ***\n" : "*** SIGABRT ***\n");

  // ============================================================================
  // Faulting address
  // ============================================================================
  // Shows WHERE the invalid memory access occurred.
  //   0x0               → null pointer dereference
  //   0x7f...           → likely heap corruption (address looks like valid heap)
  //   0x0000XXXX...     → likely stack corruption (invalid address formed by
  //                        misinterpreting non-pointer data as an address)
  WriteStr(file_desc, "Faulting address (si_addr): ");
  // printing void* as hex
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  WriteHex(file_desc, reinterpret_cast<u64>(info->si_addr));
  WriteStr(file_desc, "\n");

  // ============================================================================
  // Signal code
  // ============================================================================
  // SEGV_MAPERR (1) = address not mapped (null deref, wild pointer)
  // SEGV_ACCERR (2) = address mapped but no permission (read-only page)
  WriteStr(file_desc, "Signal code    (si_code):  ");
  WriteInt(file_desc, info->si_code);
  WriteStr(file_desc, "\n");

  // ============================================================================
  // Thread ID
  // ============================================================================
  // Match this against the crash slot thread IDs to identify the worker.
  WriteStr(file_desc, "Thread ID:                 ");
  // pthread_t is a pointer on macOS
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  WriteHex(file_desc, reinterpret_cast<u64>(pthread_self()));
  WriteStr(file_desc, "\n");

  // ============================================================================
  // Instruction pointer
  // ============================================================================
  // The exact machine instruction that caused the crash.
  // Resolve with: addr2line -e ./Lancet2 -f <address>
#if defined(__linux__) && defined(__x86_64__)
  auto* uctx = static_cast<ucontext_t*>(ctx);
  auto const rip = static_cast<u64>(uctx->uc_mcontext.gregs[REG_RIP]);
  WriteStr(file_desc, "Instruction ptr (RIP):     ");
  WriteHex(file_desc, rip);
  WriteStr(file_desc, "\n");
#elif defined(__linux__) && defined(__aarch64__)
  auto* uctx = static_cast<ucontext_t*>(ctx);
  auto const prog_ctr = static_cast<u64>(uctx->uc_mcontext.pc);
  WriteStr(file_desc, "Instruction ptr (PC):      ");
  WriteHex(file_desc, prog_ctr);
  WriteStr(file_desc, "\n");
#elif defined(__APPLE__) && defined(__x86_64__)
  auto* uctx = static_cast<ucontext_t*>(ctx);
  auto const rip = static_cast<u64>(uctx->uc_mcontext->__ss.__rip);
  WriteStr(file_desc, "Instruction ptr (RIP):     ");
  WriteHex(file_desc, rip);
  WriteStr(file_desc, "\n");
#elif defined(__APPLE__) && defined(__aarch64__)
  auto* uctx = static_cast<ucontext_t*>(ctx);
  auto const prog_ctr = static_cast<u64>(uctx->uc_mcontext->__ss.__pc);
  WriteStr(file_desc, "Instruction ptr (PC):      ");
  WriteHex(file_desc, prog_ctr);
  WriteStr(file_desc, "\n");
#else
  static_cast<void>(ctx);
  WriteStr(file_desc, "Instruction ptr:           (unavailable on this arch)\n");
#endif

  // ============================================================================
  // Backtrace
  // ============================================================================
  // Best-effort stack trace.  May produce 0 frames if the stack is
  // too corrupted to walk.  The instruction pointer above is always
  // available regardless.
  WriteStr(file_desc, "\nBacktrace:\n");
  std::array<void*, MAX_BACKTRACE_FRAMES> frames{};
  int const depth = backtrace(frames.data(), MAX_BACKTRACE_FRAMES);
  if (depth > 0) {
    backtrace_symbols_fd(frames.data(), depth, file_desc);
  } else {
    WriteStr(file_desc, "  (no frames — stack may be corrupted)\n");
  }

  // ============================================================================
  // Per-thread crash context (§2)
  // ============================================================================
  DumpAllThreadContexts(file_desc);

  // ============================================================================
  // All-thread IPs (§3, Linux only)
  // ============================================================================
#ifdef __linux__
  DumpAllThreadIPs(file_desc);
#endif

  WriteStr(file_desc, "\n============================================\n");
  WriteStr(file_desc, "Resolve addresses with:\n");
  WriteStr(file_desc, "  addr2line -e ./Lancet2 -f <address>\n");
  WriteStr(file_desc, "============================================\n");

  // Re-raise with default handler to produce a core dump.
  // raise() return value irrelevant during crash
  // NOLINTNEXTLINE(cert-err33-c)
  static_cast<void>(raise(sig));
}

}  // namespace

// ============================================================================
// §5  Public API — called from main() and AsyncWorker::Process()
//
// Lifecycle:
//   1. main()                 → InstallCrashHandler()        (once, before threads)
//   2. AsyncWorker::Process() → RegisterThreadSlot()         (once per thread)
//   3. per-window loop        → SetSlotWindowInfo()          (before ProcessWindow)
//   4. per-window loop        → ClearSlotWindowInfo()        (after ProcessWindow)
//   5. AsyncWorker::Process() → UnregisterThreadSlot()       (thread exit)
// ============================================================================

namespace lancet::base {

void InstallCrashHandler() {
  // Allocate the alternate signal stack.  This is a separate stack used ONLY
  // by the signal handler.  Without this, a stack overflow crash would corrupt
  // the very stack the handler needs to run on.
  //
  // The buffer is allocated once via a static vector (never freed, by design)
  // because the alt-stack must remain valid for the lifetime of the process
  // and SIGSTKSZ is no longer a compile-time constant on glibc >= 2.34.
  static auto alt_stack_mem = std::vector<char>(GetAltStackSize(), '\0');
  stack_t alt_stack{};
  alt_stack.ss_sp = static_cast<void*>(alt_stack_mem.data());
  alt_stack.ss_size = alt_stack_mem.size();
  alt_stack.ss_flags = 0;
  sigaltstack(&alt_stack, nullptr);

  // Register the handler for SIGSEGV and SIGABRT.
  //   SA_SIGINFO  — provides siginfo_t (faulting address) and ucontext (RIP)
  //   SA_ONSTACK  — run on the alternate stack allocated above
  //   SA_RESETHAND — auto-reset to SIG_DFL before entering the handler
  //                  (prevents infinite recursion if the handler itself crashes)
  struct sigaction act{};
  act.sa_sigaction = CrashHandler;
  act.sa_flags = SA_ONSTACK | SA_RESETHAND | SA_SIGINFO;
  sigemptyset(&act.sa_mask);
  sigaction(SIGSEGV, &act, nullptr);
  sigaction(SIGABRT, &act, nullptr);
}

auto RegisterThreadSlot() -> CrashSlotIdx {
  auto const tid = pthread_self();
  for (usize idx = 0; idx < MAX_CRASH_SLOTS; ++idx) {
    auto& slot = g_crash_slots[idx];
    // Atomic CAS: try to claim this slot (0 → 1).  Only one thread succeeds.
    int expected = 0;
    if (__atomic_compare_exchange_n(&slot.mInUse, &expected, 1, false, __ATOMIC_SEQ_CST,
                                    __ATOMIC_SEQ_CST)) {
      slot.mThreadId = tid;
      __atomic_store_n(&slot.mActive, 0, __ATOMIC_RELEASE);
      slot.mGenomeIdx = 0;
      slot.mRegion[0] = '\0';
      return idx;
    }
  }
  return INVALID_CRASH_SLOT;
}

void UnregisterThreadSlot(CrashSlotIdx const slot_idx) {
  if (slot_idx >= MAX_CRASH_SLOTS) return;
  auto& slot = g_crash_slots[slot_idx];
  __atomic_store_n(&slot.mActive, 0, __ATOMIC_RELEASE);
  __atomic_store_n(&slot.mInUse, 0, __ATOMIC_RELEASE);
}

void SetSlotWindowInfo(CrashSlotIdx const slot_idx, u64 const genome_idx, char const* region_str) {
  if (slot_idx >= MAX_CRASH_SLOTS) return;
  auto& slot = g_crash_slots[slot_idx];

  // Write data BEFORE marking active.  The release fence on the mActive store
  // ensures the signal handler sees consistent genome_idx and region data.
  slot.mGenomeIdx = genome_idx;

  usize pos = 0;
  if (region_str != nullptr) {
    for (; pos < MAX_REGION_LEN - 1 && region_str[pos] != '\0'; ++pos) {
      slot.mRegion[pos] = region_str[pos];
    }
  }
  slot.mRegion[pos] = '\0';

  __atomic_store_n(&slot.mActive, 1, __ATOMIC_RELEASE);
}

void ClearSlotWindowInfo(CrashSlotIdx const slot_idx) {
  if (slot_idx >= MAX_CRASH_SLOTS) return;
  __atomic_store_n(&g_crash_slots[slot_idx].mActive, 0, __ATOMIC_RELEASE);
}

}  // namespace lancet::base

#else  // Non-POSIX platforms (Windows, etc.) — all functions are no-ops.

namespace lancet::base {

void InstallCrashHandler() {}

auto RegisterThreadSlot() -> CrashSlotIdx {
  return INVALID_CRASH_SLOT;
}
void UnregisterThreadSlot(CrashSlotIdx /*slot_idx*/) {}
void SetSlotWindowInfo(CrashSlotIdx /*slot_idx*/, u64 /*genome_idx*/, char const* /*region_str*/) {}
void ClearSlotWindowInfo(CrashSlotIdx /*slot_idx*/) {}

}  // namespace lancet::base

#endif
