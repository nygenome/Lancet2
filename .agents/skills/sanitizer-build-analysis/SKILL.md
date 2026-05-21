---
name: sanitizer-build-analysis
description: Use when investigating a crash, intermittent test failure, suspected data race, memory error, or any sanitizer-detectable bug, AND the goal is to land a fix with a regression test. Trigger on "this crashes", "race condition", "use-after-free", "ASan says", "TSan reports", "MSan", "UBSan", "intermittent flake in tests/". Drives pixi-managed sanitizer build trees (cmake-build-asan, cmake-build-msan, cmake-build-tsan, cmake-build-ubsan), reproduces the failure, hands the trace to the sanitizer-expert agent for analysis, and lands the minimum fix with a regression test. Per-sanitizer detail (what each catches, runtime options, common findings, the mimalloc-coverage caveat) lives in references/sanitizer_matrix.md, loaded on demand. For analysis-only of an existing report, use the sanitizer-expert subagent instead.
allowed-tools: Read, Glob, Grep, Edit, Write, Bash
---

# Sanitizer build analysis on Lancet2

This skill is the end-to-end procedure for taking a sanitizer-detectable bug from "I have a reproduction" to "merged with regression test." It composes with the `sanitizer-expert` subagent for analysis; this skill is the workflow the calling session executes.

The workflow uses pixi-managed sanitizer build tasks (`pixi run configure-asan`, `build-asan`, `test-asan`, and parallel families for MSan, TSan, UBSan). Pipeline-level sanitizer runs additionally require the `LANCET_SANITIZE_BUILD` source guard in `src/lancet/main.cpp`, `benchmarks/main.cpp`, and `CMakeLists.txt` — without that guard, sanitizer builds either fail to link or run silently because Lancet2's mimalloc allocator override bypasses the sanitizer interceptors. See "Setup verification" below.

## Handbook framing

Sanitizers are **per-bug-class tools**, not "always on" tools. The right question is not "which sanitizer should I always run" — it's "which class of bug am I investigating, and which sanitizer catches that class."

### Symptom-to-sanitizer matrix (generic)

| Symptom signature | Use |
|---|---|
| Segfaults, heap corruption, weird iterator/reference crashes, "crash much later than the bug" | **ASan** |
| Optimization-only weirdness, signed overflow suspicion, casts, alignment, impossible control flow | **UBSan** |
| Nondeterministic test failures, state corruption only under load, bugs that disappear with logging, lock-free or weakly-synchronized code | **TSan** |
| Random branch choices, nondeterministic protocol fields, "garbage" values, partially-initialized structs crossing boundaries | **MSan** |
| Steadily increasing memory in tests, ownership cleanup uncertainty, exception-path leaks | **LSan** (bundled with ASan on Linux) |

Findings often come in families: ASan + UBSan together catch most memory bugs; TSan + ASan together catch races that mask memory corruption. After a finding fires under one sanitizer, re-run under adjacent sanitizers when the symptom suggests it.

### Lancet2-specific bug-class checklist

| Lancet2 pattern | Likely sanitizer |
|---|---|
| HTSlib iterator state shared across worker threads (per-thread requirement violated) | **TSan** |
| `absl::Span<u8 const>` or `std::string_view` returned over destroyed local storage | **ASan** (heap- or stack-use-after-return) |
| Iterator invalidation after `std::vector::push_back` / `reserve` in graph-traversal hot loops | **ASan** |
| Variant-store dedup hash race producing silently-wrong VCF output | **TSan** |
| `bool done` shutdown flag on worker threads without `std::atomic` or memory barrier | **TSan** |
| BAM iterator state recycled across regions without proper cleanup | **ASan** (or LSan for the leaked half) |
| Signed integer overflow in coordinate arithmetic at chromosome boundaries | **UBSan** |
| `std::shared_ptr<T>` ref-count safe but mutating `T` from multiple threads | **TSan** |
| Use-after-move: stale `data()` pointer used after the source vector was moved-from | **ASan** |

### Why sanitized builds behave differently

Instrumentation changes object layout, allocator behavior, alignment, stack frames, timing, inlining, memory reuse, and thread scheduling. Bugs that disappear under sanitizers, or appear only under sanitizers, are not noise — they mean the program depends on undefined or fragile behavior, and the tool changed the environment enough to reveal it. "ASan made it crash" often really means "ASan stopped your bug from silently corrupting memory."

## Setup verification (one-time, per repository)

Two things must be in place — both are part of the project's standard configuration:

1. **Sanitizer pixi tasks present.** `pixi task list | grep -E '(asan|msan|tsan|ubsan)'` should list 12 task names: `configure-{asan,msan,tsan,ubsan}`, `build-{asan,msan,tsan,ubsan}`, `test-{asan,msan,tsan,ubsan}`. If any are missing, the project's `pixi.toml` is out of date for sanitizer support — surface to user and stop.

2. **`LANCET_SANITIZE_BUILD` source guard present (only required for pipeline-level runs).** `grep -l LANCET_SANITIZE_BUILD CMakeLists.txt src/lancet/main.cpp benchmarks/main.cpp` — all three should match. The guard gates Lancet2's mimalloc allocator override behind a CMake option that the sanitizer pixi tasks set; without it, sanitizer builds of the `Lancet2` pipeline binary either fail to link (`multiple definition of malloc`) or run but report nothing (mimalloc bypasses the sanitizer interceptors).

If you only need to sanitize unit tests (`pixi run test-asan`, etc.), the source guard is NOT required — `TestLancet2` does not link mimalloc.

## Step 1 — Pick the right sanitizer

Match the symptom to the matrix above. The matrix in `references/sanitizer_matrix.md` covers each in detail; the short version:

- **ASan** (the default first try): crashes, heap-buffer-overflow, use-after-free, stack-buffer-overflow. Linux ASan also includes LSan. Cheap (2× CPU); always combined with UBSan in the project's `configure-asan` task.
- **TSan**: intermittent failures, hangs, output that varies between runs even with deterministic input. Roughly 5-15× CPU; needs multiple threads (use `--num-threads $(nproc)` even for normally-2-thread bugs).
- **MSan**: suspected uninitialized-memory reads in pure Lancet2 code paths. **Use sparingly** — MSan produces false positives on every call into vendored deps (htslib, abseil, spdlog, SPOA, minimap2) because those aren't MSan-instrumented.
- **UBSan-standalone**: inconsistent arithmetic, "shift exponent too large", "signed integer overflow", alignment-dependent failures. Already bundled into `configure-asan`; reach for `configure-ubsan` only when ASan's allocator interposition is interfering.

If the symptom doesn't clearly match one, start with ASan. It catches the broadest class and has the best output quality. Re-run under adjacent sanitizers if the first run is clean but the bug persists.

## Step 2 — Build under the chosen sanitizer

```bash
pixi run build-asan   # or build-msan, build-tsan, build-ubsan
```

The `build-X` task depends on `configure-X`, so a fresh tree configures and builds in one command. Subsequent runs incrementally rebuild — pixi's dependency tracker skips configure if `cmake-build-X/CMakeCache.txt` exists.

The `cmake-build-asan/`, `cmake-build-msan/`, `cmake-build-tsan/`, `cmake-build-ubsan/` directories are matched by the protected-paths hook's `cmake-build-*` glob — agent writes are blocked but build-system writes pass through. Read access from the agent into these trees is unrestricted, so reading vendored-dep source under `cmake-build-asan/_deps/` for triage works as expected.

## Step 3 — Reproduce under the sanitizer

For unit-test reproduction:

```bash
pixi run test-asan -- "[failing test name]"   # or test-msan/test-tsan/test-ubsan
```

The `--` separates the pixi task args from arguments forwarded to `TestLancet2`. The trailing `[tag]` filter limits the run to one Catch2 test, which is what you want during triage.

For pipeline-level reproduction (requires the `LANCET_SANITIZE_BUILD` source guard):

```bash
pixi run ./cmake-build-asan/Lancet2 pipeline \
  --tumor $LANCET_TEST_SOMATIC_TUMOR \
  --normal $LANCET_TEST_SOMATIC_NORMAL \
  --reference $LANCET_TEST_SOMATIC_REFERENCE \
  --region $LANCET_TEST_SOMATIC_REGION_SMALL \
  --num-threads $(nproc) \
  --out-vcfgz /tmp/triage.vcf.gz \
  2> /tmp/sanitizer.log
```

The `pixi run` prefix is required for sanitizer pipeline runs because sanitizer builds are dynamically linked (`-DLANCET_BUILD_STATIC=OFF`); the sanitizer runtime libraries (`libasan.so`, `libtsan.so`, etc.) and `libc++` resolve from the pixi env's library path. This differs from production and profiling builds, which are statically linked.

Use `_REGION_SMALL` (1 Mb) regions to keep iteration fast — sanitizer-instrumented binaries are 2-15× slower depending on which sanitizer. Save the sanitizer output to a file (`2> /tmp/sanitizer.log`) — output is verbose and easy to lose to terminal scrollback.

If the bug doesn't reproduce on first try:
- TSan: increase thread count, run multiple times (10×) to catch intermittents.
- ASan/UBSan: try `ASAN_OPTIONS=detect_stack_use_after_return=1:check_initialization_order=1` for additional checks.
- Any sanitizer: try a different region, or a different fixture (germline NA12878 vs somatic HCC1395).

## Step 4 — Hand off to sanitizer-expert agent

Hand the captured sanitizer output to the agent:

```
Use the sanitizer-expert subagent to analyze this output:
[paste the sanitizer report]

The reproduction is: [your command]
```

The agent reads the cited frames in `src/lancet/`, identifies the root cause, and proposes the minimum fix. It has its own memory file at `.Codex/agent-memory/sanitizer-expert.md` recording known-benign warnings and patterns of new warnings that turned out to be real bugs — consult that before reasoning from scratch.

For unfamiliar diagnostic categories ("what does this UBSan report mean?"), `references/sanitizer_matrix.md` has per-sanitizer notes including option semantics and Lancet2-specific gotchas.

## Step 5 — Write a regression test FIRST

Before implementing the fix, add a Catch2 test in `tests/<layer>/` that reproduces the bug deterministically. Use the `add-cpp-test` skill for the canonical procedure. **The test must fail under the sanitizer build before the fix is applied.**

For race conditions, deterministic reproduction is hard but not impossible. Patterns that work: constraining thread count to 2-4 (small enough for the test to be fast, large enough for contention), injecting `std::this_thread::yield()` calls at the suspected race window, using `absl::Mutex::AssertHeld` to verify locking invariants.

## Step 6 — Apply the minimum fix

In `src/lancet/<layer>/`, make the smallest change the agent proposed. The fix must address the root cause, not the symptom. Common Lancet2 fix patterns:

- Per-thread instances of a shared resource (most HTSlib iterator races).
- Explicit synchronization on a previously-implicit shared variable (most variant-store races).
- Bounds-checking on a previously-trusted index (most ASan findings).
- Proper RAII for an HTSlib resource (most leak findings).

## Step 7 — Verify under sanitizer and under release

Run the regression test under the sanitizer build to confirm the fix:

```bash
pixi run test-asan -- "[regression test name]"
```

Then verify the fix doesn't regress the release build:

```bash
pixi run build-release
pixi run test
pixi run lint-check
```

A fix that passes the sanitizer but fails clang-tidy needs further work.

## Step 8 — Suppression philosophy

If a sanitizer finding cannot be fixed (vendored-dep behavior, intentional pattern), suppression is a last resort. Use ignorelist files (`sanitizer_ignorelist.txt`) or `__attribute__((no_sanitize("...")))` deliberately, locally, with a comment explaining why the tool is wrong for this region. **Never suppress because the report is annoying.** The failure mode to avoid is "we suppressed it because triage was inconvenient."

## Step 9 — Commit

Use the conventional-commit prefix `fix:` and reference the sanitizer in the body:

```
fix: per-thread BAM iterator in read_collector

Reproduced under TSan (cmake-build-tsan, --num-threads 8) on
chr1:1000000-2000000. Race was on shared bam1_t* between worker
threads — HTSlib iterators are not thread-safe. Fix: per-worker
iterator instance with RAII cleanup.

Regression test: tests/hts/read_collector_test.cpp
```

Lancet2's `.chglog/config.yml` does not support scopes; mention the layer in the subject text instead.

## When NOT to use this skill

Do not use this skill for bugs that aren't sanitizer-detectable (logic errors, wrong VCF schema, race-free correctness gaps); use the standard `add-cpp-test` workflow. Do not use it as a routine pre-merge check — the project's CI already runs ASan on the test suite. Use this skill specifically when you have a symptom and need to find the root cause.

## When the bug is in a vendored dependency

If the agent identifies the root cause inside HTSlib, minimap2, SPOA, longdust, or another vendored library under `cmake-build-*/_deps/`, the fix path changes. Do NOT patch the vendored code (write-blocked, lost on next reconfigure). Instead, work around it at the Lancet2 callsite. Common workarounds:
- **Per-thread instances** (for thread-unsafe libraries — most HTSlib uses).
- **Defensive copying** (for libraries with unstable iterator semantics).
- **Pre-validation** (for libraries with insufficient input checking).

Document the workaround with a comment citing the upstream behavior.

## A note on mimalloc coverage

Sanitizer builds use libc malloc, not Lancet2's production mimalloc (the `LANCET_SANITIZE_BUILD` source guard makes this unconditional in sanitizer trees). A clean ASan run does NOT prove production-with-mimalloc is clean. Bugs that depend on mimalloc's allocation pattern, thread-local heap layout, or page-recycling timing will not surface here. The trade-off is intentional — mimalloc's preprocessor-level override of `malloc/free/new/delete` bypasses sanitizer interceptors entirely, so the choice is "sanitize without mimalloc" or "don't sanitize at all." The guard's rationale and exact scope are documented inline in `src/lancet/main.cpp` and `CMakeLists.txt`; read them when adjusting the sanitizer build configuration.

## References

- `references/sanitizer_matrix.md` — Per-sanitizer detail (ASan/MSan/TSan/UBSan/LSan): what each catches, build flag rationale, runtime options, common findings, Lancet2-specific notes, canonical Clang documentation URLs. Load at Step 1 (picking the sanitizer) and Step 4 (interpreting unfamiliar output).
