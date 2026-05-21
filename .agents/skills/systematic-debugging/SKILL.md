---
name: systematic-debugging
description: Use when investigating a bug, test failure, unexpected output, or unexpected behavior in Lancet2 source — before proposing or attempting any fix. Trigger on "this test is failing", "the variant counts don't match expectations", "I'm getting wrong output from the pipeline", "X stopped working after Y change", "we have a flaky test", "the same test fails for different reasons each run", "I tried a fix but it didn't work", "should I just patch around this for now", "the assertion fires but I don't understand why". Walks the four-phase root-cause-first methodology (root-cause investigation, pattern analysis, hypothesis and testing, implementation), applies the rule that no fix lands without an articulated root cause, and stops the session on the architecture-question pattern (3+ fixes failed, each revealing a new symptom in a different place). Defers per-bug-class detail: sanitizer findings to `sanitizer-build-analysis` skill and `sanitizer-expert` agent; performance to `profile-and-optimize`; probe-tracking attribution to `probe-tracking` skill and `probe-interpreter` agent; VCF schema to `vcf-validator`.
allowed-tools: Read, Glob, Grep, Bash, Edit
---

# Systematic debugging on Lancet2

The first answer to "this is broken" is not "let me try a fix and see." It is "what is the root cause?" Symptom-driven fixes mask underlying issues, accumulate technical debt, and tend to surface the same bug in a different shape weeks later. This skill walks the procedure for finding the root cause first and fixing it once.

The rule: no fix lands without an articulated root cause. If you cannot state in one sentence what triggers the bug and why, the fix you write addresses a symptom, and symptom fixes are how this codebase accumulates tech debt fastest.

The procedure has four phases. The transitions between phases are gates: complete each phase before proceeding to the next.

## Phase 1 — Root-cause investigation

Read the failure carefully before reasoning about it. The most common error in this phase is jumping past the actual evidence to a hypothesis.

**Read the error fully.** Stack traces, sanitizer reports, `LANCET_ASSERT` failures, and Catch2 SECTION names all carry diagnostic information that is often dismissed as boilerplate. Read them. Note the file:line, the variable values, the failing assertion text. The phrase "the assertion fires but I don\'t know why" is itself a signal that the assertion text was skimmed rather than read.

**Reproduce consistently.** Can you trigger the failure reliably? What command, what input, what test? If reproduction is intermittent — same input, different result — say so explicitly. Intermittency usually points to threading, memory, or environment issues (TSan / ASan / sanitizer differences) rather than logic bugs, and the appropriate next step is `sanitizer-build-analysis`, not more guessing.

**Check what changed.** `git diff HEAD~1`, `git log --since="1 week ago"`, recent-PR list. Regressions often track to a specific commit, and bisection is faster than reasoning when the regression window is narrow.

**Trace data flow when the failure is deep in the call stack.** When a wrong value surfaces in `cli/` but the real bug is in `cbdg/`, follow the value backward: where does it first appear? Who passed it forward unchecked? See `root-cause-tracing.md` in this directory for the complete backward-tracing procedure.

**Instrument boundaries when symptoms cross components.** Lancet2 has natural component boundaries: layer transitions (`base/` → `hts/` → `cbdg/` → `caller/` → `core/` → `cli/`), per-thread workers feeding the shared `VariantStore`, the producer/consumer queue between `pipeline_executor` and `async_worker`. When a symptom appears in one layer but the trigger is upstream, add `LOG_DEBUG` calls (defined in `src/lancet/base/logging.h`) at each boundary to log what enters and what exits. The Release build compiles `LOG_DEBUG` out automatically, so the instrumentation is free at production cost; iterate against `cmake-build-debug/` where the calls stay live.

A common Phase-1 anti-pattern is "the bug is probably in X" without evidence. Without evidence, the diagnosis is a guess; the fix that follows lands somewhere in X and either fixes the symptom (root cause untouched) or makes things worse (new bug introduced).

## Phase 2 — Pattern analysis

Before forming a hypothesis, find the related shape in the codebase.

**Find a working analogue.** Lancet2 has a lot of similar-shape code: graph traversals in `cbdg/`, allele assignment in `caller/`, BAM read collection in `hts/`, sharded sinks in `core/`. If your code is broken, find a piece of code that works the same way and compare. The differences pinpoint the bug.

**Read reference implementations completely.** When implementing a pattern (Welford\'s algorithm for online stats, Dirichlet-Multinomial genotype likelihoods, BCALM 2 sign continuity for k-mer canonicalization, the MaxFlow walk-tree arena, SPOA convex scoring with int16 SIMD lanes), read the cited reference end-to-end. Skimming the section that "looks relevant" is how subtle errors enter — the assumption you missed two pages above is exactly the one that does not hold in your case.

**List every difference between working and broken.** Do not dismiss differences as "that can\'t matter" — that is bias talking. Write them down. The bug is usually in a difference initially classified as insignificant.

**Understand all dependencies.** What does this code need from upstream? What invariants does it assume about its inputs? What does it produce that downstream code relies on? In Lancet2, frequently-missed invariants include: per-thread `Extractor`/`Iterator` lifetime in `hts/`, k-mer canonicalization rules in `cbdg/`, allele-score sign and zero-handling in `caller/`, and the shard-key computation in `core/VariantStore`.

## Phase 3 — Hypothesis and testing

Form a hypothesis explicitly, test minimally, verify before continuing.

**State the hypothesis in one sentence.** "I think `<X>` is the root cause because `<Y>`." Specific, testable, falsifiable. "I think the SPOA scoring is wrong" is not a hypothesis. "I think `MSA_MATCH_SCORE = 0` causes the alignment to prefer mismatches when both have score 0, which is why the variant is split into two raw_variants" is a hypothesis.

**Test minimally.** The smallest change that would falsify or confirm the hypothesis. One variable at a time. Multiple changes at once are diagnostic poison: when something breaks, you cannot tell which change broke it; when something passes, you cannot tell which change made it pass.

**Verify the hypothesis was correct.** Did the fix produce the predicted behavior? If yes, proceed to Phase 4. If no, the hypothesis was wrong; return to Phase 1 with the new evidence. Do not stack a second fix on top of an unconfirmed first one.

**When you don\'t understand something, say so.** Lancet2\'s deep parts (`cbdg/`, `caller/`, the async pipeline in `core/`) reward humility. "I don\'t understand why this lock does not deadlock" is a useful starting point; "this probably works, let me move on" is not.

## Phase 4 — Implementation

Fix the root cause. Prove the fix with a regression test.

**Write the failing test first.** A regression test that fails today and passes after the fix is the proof that the fix is real. Use the `add-cpp-test` skill for the Catch2 fixture; the test layout (`tests/<layer>/<file>_test.cpp`), the assertions, the test-binary invocation conventions, and the project\'s Catch2 idioms cheat sheet are encoded there.

**Implement the minimum fix.** Address the root cause; do not bundle "while I\'m here" cleanups. Cleanups belong in a separate commit so a future bisect can isolate the fix from cosmetic changes.

**Verify with `/fix-and-validate`.** Lint, IWYU, and the test suite must all pass against the Release tree. The PreToolUse `pre_commit_gate.sh` hook will block any commit without a fresh validation marker (a `git stash create` hash matching the about-to-be-committed working tree), so this verification step is enforced regardless.

**If the fix does not resolve the issue, stop.** Do not stack a second fix on top. Return to Phase 1 with the new information.

## When 3+ fixes have failed: question the architecture

When three different fixes targeting three different root-cause hypotheses have all failed, the diagnosis is wrong at a level deeper than any individual hypothesis. The pattern looks like:

- Each fix reveals a new symptom in a different file or layer.
- Each fix requires a wider refactor than the last.
- Fixes that should be independent interact unexpectedly when stacked.

When that pattern holds, the architecture itself may be wrong for the workload. Stop fixing. Surface the pattern to the user. Discuss whether the right move is to refactor the structure rather than continue patching symptoms. Do not attempt a fourth fix without that discussion.

This is not a hypothesis failure — it is a diagnosis failure at a higher level. The signal is strong enough that ignoring it produces worse outcomes than stopping and reconsidering.

## Common rationalizations to refuse

- "Quick fix for now, investigate later." The first fix sets the pattern; later investigation rarely happens, and the symptom-fix becomes load-bearing.
- "Time pressure means no time for process." Systematic debugging is faster than guess-and-check thrashing — every iteration of guess-and-check costs the time the investigation would have cost.
- "I\'ll write the test after confirming the fix works." Untested fixes do not stick. The regression test is the proof; without it, "the fix works" is an opinion.
- "I see the symptom, let me fix it." Seeing a symptom is not understanding the cause. Treating one is treating the other only by accident.
- "One more fix attempt." After two failed attempts targeting different hypotheses, the third without re-investigation lands as #4 in the architecture-question pattern.

## When NOT to use this skill

- **Vendored deps** under `cmake-build-*/_deps/` are read-only by `block_protected_paths.py`; the only fix path is at the Lancet2 callsite. The bug class is "how does Lancet2 use this dependency", not "what does the dependency do".
- **Configuration drift, environment issues, CI infrastructure failures.** These typically do not need root-cause source debugging — check the env, read the CI logs, fix the config or the workflow.
- **Style or linting violations.** Use `clang-tidy-discipline` for the procedural how-to (diagnose, fix or scope-suppress, never `--fix`).
- **Performance regressions where correctness is intact.** Use `profile-and-optimize` instead; the methodology is different (gperftools/pprof + benchmark suite, not assertion-trace).

## Related Lancet2 surfaces

- `add-cpp-test` skill — for writing the failing Catch2 regression test (Phase 4).
- `/fix-and-validate` slash command — runs the project\'s lint + test gates after the fix is in (Phase 4).
- `sanitizer-build-analysis` skill and `sanitizer-expert` agent — when the bug surfaces under ASan / MSan / TSan / UBSan, or when intermittency points at memory or threading.
- `assembly-and-calling-expert` agent — for `cbdg/` + `caller/` correctness reasoning when Phase 2 needs deep domain knowledge.
- `probe-tracking` skill and `probe-interpreter` agent — when the bug is "Lancet2 misses a variant"; produces a per-stage attribution from the probe-tracking pipeline so Phase 1 starts with structured evidence rather than a vague missing-variant report.
- `vcf-validator` agent — when the bug touches VCF schema (FORMAT/INFO/FILTER fields).
- `perf-analyst` agent — when the bug is performance, not correctness.
- `fresh-reviewer` agent — invoke after a non-trivial fix lands; the cold review catches what the writer rationalized away.

## Reference files in this directory

- `root-cause-tracing.md` — the backward-tracing procedure from symptom to source for failures deep in the call stack.
- `defense-in-depth.md` — adding validation at multiple tiers after the root cause is fixed, so the same bug class becomes structurally harder to reintroduce.

