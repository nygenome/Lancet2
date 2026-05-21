# Sanitizer matrix

Per-sanitizer reference for the four sanitizers Lancet2's pixi tasks
support. Each section names: what the sanitizer catches, the pixi
tasks that drive it, common findings and how to read them, and
Lancet2-specific notes (the mimalloc situation, vendored-dep coverage
gaps, etc.).

The Clang documentation is the authoritative source for sanitizer
behavior. Each section ends with the canonical doc URL — read it
before customizing options or interpreting unfamiliar output.

This reference is loaded on demand by the `sanitizer-build-analysis` skill
when configuring a build, interpreting a report, or deciding which
sanitizer fits a given symptom. Read the section that matches the
sanitizer you're using, not the whole file.

## AddressSanitizer (ASan) — the default

**What it catches.** Memory errors at runtime: heap-buffer-overflow,
heap-use-after-free, stack-buffer-overflow, stack-use-after-return
(off by default — see `ASAN_OPTIONS` below), stack-use-after-scope,
global-buffer-overflow, use-after-poison, double-free, invalid-free,
container-overflow (with annotated containers), initialization-order-
fiasco. On Linux, also catches leaks via the bundled LeakSanitizer
(LSan) — controlled via the same `ASAN_OPTIONS` namespace.

**Pixi tasks.** `pixi run configure-asan`, `pixi run build-asan`,
`pixi run test-asan`. The configure includes UBSan in the same tree
because the combination is cheap and catches more bugs than ASan
alone.

**When to reach for it.** This is the default first sanitizer. Crash
output with `SEGV`, `heap-buffer-overflow`, `use-after-free`, or
`stack-buffer-overflow` in the symptom is ASan territory. Slow
accumulation of memory over a long run (suspected leak) is also ASan,
via its bundled LSan.

**Overhead.** Roughly 2× CPU, 2-3× memory. Acceptable for
`TestLancet2` runs and for `LANCET_TEST_*_REGION_SMALL` pipeline runs;
the `_SMALL` regions are 1 Mb and finish in seconds even instrumented.
Whole-chr regions take long enough that you'll want to reduce thread
count or pick a smaller region.

**Most relevant runtime options** (set via `ASAN_OPTIONS`,
colon-separated):

| Option | Effect |
|:-------|:-------|
| `halt_on_error=1` | Abort on first error rather than continue (default 1; turn off only when collecting multiple errors) |
| `detect_leaks=1` | Enable LSan integration on Linux (default 1 on Linux) |
| `detect_stack_use_after_return=1` | Catch UAR (off by default; significant overhead, but the failures are nasty) |
| `abort_on_error=1` | `abort()` instead of `_exit()` so debuggers and core-dump tooling see the failure |
| `print_stacktrace=1` | Always print stack on error (default 1) |
| `symbolize=1` | Resolve symbols (default 1; if traces show raw addresses, `llvm-symbolizer` is missing from PATH) |
| `strict_string_checks=1` | More aggressive string-function bounds checking (catches more, more FPs) |
| `check_initialization_order=1` | Detect initialization-order fiasco between TUs |
| `log_path=/tmp/asan` | Write reports to `/tmp/asan.<pid>` instead of stderr (useful for long runs) |
| `disable_coredump=0` | Allow core dumps on abort (off by default for binary-size reasons) |

A useful starting set for triage: `ASAN_OPTIONS=halt_on_error=1:abort_on_error=1:detect_stack_use_after_return=1:check_initialization_order=1:strict_string_checks=1`.

**Lancet2-specific notes.**

- Sanitizer builds use libc malloc, NOT mimalloc (the
  `LANCET_SANITIZE_BUILD` macro disables the mimalloc include and
  link in sanitizer trees; the macro's rationale is documented
  inline in `src/lancet/main.cpp` and `CMakeLists.txt`). A clean
  ASan run does NOT prove production-with-mimalloc is clean. Bugs
  that depend on mimalloc's allocation pattern, thread-local heap
  layout, or page recycling will not surface here.
- Vendored deps under `cmake-build-*/_deps/` ARE built with
  ASan instrumentation when present in the same tree (htslib,
  abseil, spdlog, SPOA, minimap2, longdust). Findings inside these
  deps are real ASan findings — not false positives — but the fix
  path is "work around at the Lancet2 callsite," not "patch the dep."
  See the SKILL's "When the bug is in a vendored dependency" section.
- LSan integration is on by default on Linux. If a leak surfaces,
  the report includes the allocation stack; investigate the call site
  in `src/lancet/` first — most leaks are RAII gaps in HTSlib
  resource handling.

**Clang doc.** https://clang.llvm.org/docs/AddressSanitizer.html


## MemorySanitizer (MSan) — uninit-mem detection (use sparingly)

**What it catches.** Reads of uninitialized memory across function
boundaries: stack variables read before assignment, heap memory read
before initialization, structure padding read as data, and so on.
The tracking is byte-level: an allocation starts entirely
uninitialized, becomes initialized as you write to it, and reading
any uninitialized byte is a finding.

**Pixi tasks.** `pixi run configure-msan`, `pixi run build-msan`,
`pixi run test-msan`. The configure includes
`-fsanitize-memory-track-origins=2` so the report tells you not just
that a byte was read uninitialized but also where it was allocated.
This is the slow form (track-origins=2 is roughly 2× the overhead of
plain MSan) but the practical form for triage.

**When to reach for it.** Symptoms strongly suggesting uninitialized
reads: results that vary across runs with the same input, stack-
variable values that sometimes-zero-sometimes-garbage, conditional
branches that depend on bits that should never have been set. ASan
will catch some uninitialized reads if the memory is also out of
bounds, but for in-bounds-but-uninitialized, ASan is silent. MSan is
the right tool.

**Overhead.** Roughly 3× CPU baseline, 6× with `track-origins=2`.
Memory overhead is comparable to ASan. Use the `_SMALL` test regions
exclusively.

**Hard caveat — vendored deps and false positives.** MSan requires
EVERY linked library to be MSan-instrumented, including libc++.
Lancet2's vendored deps (htslib, abseil, spdlog, SPOA, minimap2,
longdust) are built without MSan instrumentation by the standard
configure-msan task because their build systems don't support the
flag straightforwardly. As a result, MSan considers every byte
returned from those libraries uninitialized, and you get a flood of
false positives on any code path that touches them.

In practice this means MSan is only useful for debugging suspected
uninit-mem in pure Lancet2 code paths that don't call into vendored
deps. For the typical Lancet2 hot path (graph construction calls
into htslib and minimap2 constantly), MSan is unusable as-is.
ASan should be the default; reach for MSan only when ASan and TSan
have come up empty AND the symptom strongly suggests uninit-mem
AND the suspect code is in a pure-Lancet2 region.

**If you genuinely need MSan on a code path that touches vendored
deps** the only honest path is rebuilding those deps with MSan
instrumentation, which is not a small undertaking. Document the
intent in `notes/<feature>/msan_setup.md` before starting.

**Most relevant runtime options** (`MSAN_OPTIONS`):

| Option | Effect |
|:-------|:-------|
| `halt_on_error=1` | Abort on first error |
| `print_stats=1` | Print summary stats at exit |
| `poison_in_dtor=1` | Mark destroyed objects as poisoned for use-after-destroy detection |
| `wrap_signals=1` | Intercept signal handlers (default 1) |
| `report_umrs=1` | Report use-of-uninitialized values (default 1) |

**Clang doc.** https://clang.llvm.org/docs/MemorySanitizer.html


## ThreadSanitizer (TSan) — race detection

**What it catches.** Data races (one write + one access on the same
memory location, neither protected by mutual exclusion or atomics),
mutex misuse (double-lock, unlock-without-lock, signal handlers
calling unsafe code), atomicity violations, and a class of
"happens-before" violations. TSan is precise: a finding is almost
always a real race, though sometimes a benign one.

**Pixi tasks.** `pixi run configure-tsan`, `pixi run build-tsan`,
`pixi run test-tsan`.

**When to reach for it.** Intermittent test failures, flakiness in
the CI pipeline that disappears on re-run, hangs that resolve
themselves, output that varies between runs even with deterministic
inputs. Lancet2 is heavily multi-threaded (window-batch parallelism,
per-thread genotyper state, shared variant store) and races are
the second-most-common bug class after memory errors.

**Overhead.** Roughly 5-15× CPU and 5-10× memory — much higher than
ASan. Use the `_SMALL` test regions exclusively for pipeline-level
runs. For unit tests, raw TSan is usually fine. Importantly,
**increase thread count** when reproducing under TSan: races often
require contention to surface, and a single-threaded TSan run sees
nothing. Use `--num-threads $(nproc)` even if the bug originally
reported with 2 threads.

**Most relevant runtime options** (`TSAN_OPTIONS`):

| Option | Effect |
|:-------|:-------|
| `halt_on_error=1` | Abort on first race |
| `second_deadlock_stack=1` | Print stacks for both lock acquisitions in deadlock reports (essential) |
| `history_size=7` | Bigger access history → catches more races at the cost of memory (default 2; 7 is the max) |
| `report_destroy_locked=1` | Report mutex destroyed while held (default 1) |
| `report_thread_leaks=1` | Report leaked threads (default 1) |
| `flush_memory_ms=N` | Flush memory state every N ms; use for long runs to keep memory bounded |

A useful starting set: `TSAN_OPTIONS=halt_on_error=1:second_deadlock_stack=1:history_size=7`.

**Lancet2-specific notes.**

- Per-thread BAM iterator pattern (see `src/lancet/hts/iterator.h`
  and `src/lancet/hts/extractor.h`, plus the per-sample
  `SampleExtractors` map in `src/lancet/core/read_collector.h`):
  HTSlib's iterator state is NOT thread-safe. Lancet2 owns one
  `Iterator` per worker, fed from the per-thread `Extractor`. TSan
  findings here usually mean the per-thread invariant has been
  broken — a shared state path was added that bypasses the
  per-thread instance.
- Variant-store sharding (`src/lancet/core/variant_store.h`,
  256 shards): the sharding is the synchronization primitive. TSan
  findings on variant insertion are usually a wrong shard-index
  computation, not a missing lock.
- The pipeline's outer windowing loop is NOT a race source — it
  uses a `concurrentqueue` MPMC FIFO that's lock-free. Findings
  there are typically inside the per-thread genotyper state, which
  should never be shared.
- TSan does NOT understand `mimalloc`'s thread-local heap (which
  doesn't matter in sanitizer builds since we use libc malloc) but
  it also doesn't understand `absl::Mutex`'s internal state. Findings
  inside `absl::Mutex` itself are almost always false positives;
  findings on the protected data are real.

**Clang doc.** https://clang.llvm.org/docs/ThreadSanitizer.html


## UndefinedBehaviorSanitizer (UBSan)

**What it catches.** Specific categories of undefined behavior:
signed-integer overflow, unsigned-integer overflow (under
`-fsanitize=integer`, opt-in), shift-exponent-too-large, alignment
violations, divide-by-zero, null-pointer dereference (`-fsanitize=null`),
type-mismatch (downcasting to wrong type), invalid `bool` value,
invalid enum value (under `-fsanitize=enum`), array bounds (with
`-fsanitize=bounds`), and a dozen others. Each check is independently
toggleable.

**Pixi tasks.** Two ways to use UBSan:
- Combined with ASan (the default in `pixi run configure-asan`):
  `-fsanitize=address,undefined`. Cheap; you get both for the price
  of ASan.
- Standalone (`pixi run configure-ubsan`, `build-ubsan`, `test-ubsan`):
  `-fsanitize=undefined` only. Much faster than ASan-combined; useful
  when ASan's allocator interposition is interfering with the bug
  you're chasing, or when you want a quick UBSan-only sweep over the
  test suite.

**When to reach for standalone UBSan.** Symptoms involving
inconsistent arithmetic, unexpected sign changes, "the same input
produces different output on different machines" (alignment), or
diagnostics that mention "shift exponent" or "signed integer overflow."
For most work, the ASan-combined form is fine.

**Overhead.** Negligible to moderate (~10-30%) depending on which
checks are enabled.

**Most relevant runtime options** (`UBSAN_OPTIONS`):

| Option | Effect |
|:-------|:-------|
| `halt_on_error=1` | Abort on first UB |
| `print_stacktrace=1` | Print stack on each finding (off by default — turn on for triage) |
| `report_error_type=1` | Include the UB category name in the report (default 1) |
| `silence_unsigned_overflow=1` | Suppress unsigned-overflow reports (NOT UB, but `-fsanitize=integer` reports it; usually noise) |

A useful starting set: `UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1`.

**Lancet2-specific notes.**

- Lancet2 uses signed integers for read positions, k-mer indexes, and
  graph node IDs. UBSan's signed-overflow check is genuinely useful
  here — overflow at negative signed values has bitten the codebase
  before (see fresh-reviewer's agent-memory if populated).
- Alignment findings on hot paths are usually misaligned reads from
  packed BAM records. The fix is `memcpy` to a stack temporary, not
  pointer cast.
- The SIMD code paths in `src/lancet/base/repeat.cpp`
  (`IsWithinHammingDist`, AVX2 + NEON intrinsics) do bit-level work
  that some UBSan checks may flag spuriously. If a finding lands
  inside the intrinsic-laden body of one of those functions and the
  surrounding logic is correct, suppress narrowly with
  `__attribute__((__no_sanitize__("undefined")))` on that specific
  function and document the reason inline. The attribute is not
  used anywhere in the current source; if you reach for it, you are
  the first to need it and your inline comment will set the
  precedent.

**Clang doc.** https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html


## LeakSanitizer (LSan) — bundled with ASan on Linux

**What it catches.** Memory leaks at process exit: heap allocations
with no live pointer reaching them. Reports the allocation stack so
you can find the missing free.

**How to use it.** LSan is bundled with ASan on Linux — no separate
build. Set `ASAN_OPTIONS=detect_leaks=1` (the default on Linux). For
standalone leak-only mode (no other ASan checks), use
`-fsanitize=leak` instead of `-fsanitize=address` — much cheaper.
Lancet2's pixi tasks don't ship a leak-only configure; the ASan tree
covers it.

**When to reach for it.** Slow memory growth over a long pipeline run.
Run with `ASAN_OPTIONS=detect_leaks=1` and let the pipeline finish
or interrupt with `SIGINT` (LSan reports at exit either way).

**Lancet2-specific notes.**

- Most leaks are HTSlib resource gaps: `bam1_t`, `htsFile*`,
  `bcf_hdr_t*`, `tbx_t*`. These have project-internal RAII wrappers in
  `src/lancet/hts/` — see the `hts.md` rule for the lifetime contract.
- Apparent "leaks" on shutdown that come from spdlog or absl static
  registries are intentional and not findings — these libraries
  deliberately leak singletons rather than racing against destruction.
  LSan's default suppressions list handles most of these; if a stack
  in your report ends in `__cxa_atexit` or similar, ignore it.

**Clang doc.** https://clang.llvm.org/docs/LeakSanitizer.html
