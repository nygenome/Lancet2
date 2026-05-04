---
name: add-cpp-test
description: Use when adding a unit or integration test, writing a new Catch2 TEST_CASE, scaffolding a test fixture, starting TDD work, OR figuring out how to run / filter / debug Catch2 tests flexibly. Trigger on "add a test", "write a test for", "I need test coverage", "TDD this", "run just this test", "how do I filter tests by tag", "Catch2 generators", "Catch2 sections", "Catch2 matchers", "reproduce a random-order failure". Enforces Catch2 v3.14.0 patterns, the one-test-one-implementation TDD discipline, the layer-direction rule for test sources, and the tests/<layer>/<file>_test.cpp layout convention. The Catch2 idiom cheat sheet and the bench-binary invocation reference live as separate files under `references/`, loaded on demand when the workflow steps point to them.
allowed-tools: Read, Glob, Grep, Edit, Write, Bash
---

# Adding a test to Lancet2

This skill formalizes the canonical procedure for adding a Catch2 test. The discipline matters because Claude's default is implementation-first; under context pressure it may rewrite a test to match buggy code rather than fix the code. The one-test-one-implementation rhythm prevents that failure mode.

## Step 1 — Locate the right test file

Tests live in `tests/`. Mirror the structure of `src/lancet/`: a test for `src/lancet/caller/genotyper.cpp` lives at `tests/caller/genotyper_test.cpp` (the convention is the `_test.cpp` suffix; verify by globbing `tests/`). If the unit under test does not yet have a test file, create one matching the existing naming convention. If it does, append to the existing file.

## Step 2 — Read existing tests in the same area

Before writing, read at least two existing tests in the relevant subdirectory in full. Catch2 idioms used across the project will reveal themselves: the choice between `TEST_CASE` with `SECTION` versus `TEMPLATE_TEST_CASE`, the convention for fixture setup, the use of `REQUIRE` versus `CHECK`, and the project's tolerance conventions for floating-point comparison.

## Step 3 — Write ONE failing test

Write a single `TEST_CASE` for the specific behavior you intend to exercise. Do not write multiple tests speculatively; under context pressure you will rewrite them to match the implementation rather than reflect the actual specification. The test should fail for the expected reason — a `REQUIRE` that compares against the expected output, with a clear failure message.

For correctness-critical math (genotype likelihoods, Phred conversions, scoring), include numerical-equivalence checks against a known-correct reference value computed offline. Use `Approx` (or the project's tolerance helper) with explicit absolute and relative tolerances.

## Step 4 — Build and run, confirm failure

Run `pixi run test`. The output should show your new test failing for the expected reason. If it fails for a different reason (compilation error, wrong fixture, unrelated bug), fix that first; you cannot trust an implementation against a test that does not fail correctly.

## Step 5 — Implement the minimum change

In `src/lancet/<layer>/`, make the smallest change that turns the test green. Respect the layer-direction rule (the at-write hook will catch violations, but knowing the rule prevents wasted effort). Do not refactor surrounding code; resist the temptation to add helpers, extract methods, or improve names. Keep this commit narrow.

## Step 6 — Build and run, confirm green

Run `pixi run test` again. The new test must pass. All previously-passing tests must continue to pass. If any other test fails, the change has unintended scope; revert and reconsider.

## Step 7 — Format and lint

The PostToolUse hook should have already run `pixi run fmt-fix`. Run `pixi run lint-check` manually before commit to catch any clang-tidy violations the at-write naming hook missed.

## Step 8 — Commit

Use the project's chglog-filtered commit prefix. Test-only commits use `chore:` (the catch-all that renders as "Refactoring" in CHANGELOG.md, since the project's `.chglog/config.yml` does not filter on `test`). A test that lands together with the code it exercises uses `fix:` or `feat:` matching the source change. Keep the commit message imperative and specific: `chore: add Phred-overflow test for posterior_base_qual` is better than `add tests`.

## Step 9 — Stop here

Resist the temptation to add more tests in the same change. Each additional test is a separate commit and ideally a separate session. The vertical-slice discipline is what makes the workflow reliable; widening scope is how it breaks down.

## When NOT to use this skill

Do not use this skill for benchmarks; benchmarks live in `benchmarks/` and the procedure is different (use the `profile-and-optimize` skill, which incorporates the canonical benchmark workflow). Do not use it for end-to-end pipeline tests that exercise the binary against real BAMs; that is the `/e2e-pipeline-test` slash command's job.

## Layer-direction reminder

Test sources in `tests/<layer>/` may include from any layer of `src/lancet/`. Tests are exempt from the layer-direction rule because they exist to exercise the layered code, not to be part of it.

## References

These are loaded on demand when this skill points to them.

- `references/catch2_idioms.md` — Catch2 v3.14.0 cheat sheet: assertion families, `TEST_CASE`/`SECTION`, tags, type-parametrised tests, generators, fixtures, matchers, logging, runtime skipping, plus the note about why benchmarks are out of scope. Read this when writing or extending a test (Steps 2 and 3).
- `references/catch2_invocation.md` — Running the test binary directly: filtering by name/tag, running a specific section or generator value, discovery (list tests/tags/reporters), reproducing a failed random-order run, debugging modes (break/abort/no-throw), reporters, sharding, and a recipe table for common scenarios. Read this when iterating on a failing test (Steps 4 and 6).
