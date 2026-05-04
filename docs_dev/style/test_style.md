# Lancet2 test style

This document is the project-wide reference for adding or editing tests
under `tests/`. It is a companion to `cpp_style.md` (code style and
clang-tidy conventions) and `code_comments.md` (comment conventions).
Where rules in those documents apply equally to source and test code,
they are not repeated here.

The rules below apply to every layer (`base`, `hts`, `cbdg`, `caller`,
`core`, `cli`) and to every kind of test (unit, integration, property,
cross-validation, golden-master). Microbenchmarks live in
`benchmarks/` (google/benchmark target) and are out of scope for this
guide — see the `profile-and-optimize` skill for benchmark conventions.

---

## 1. File layout

- One test file per source unit: `tests/<layer>/<unit>_test.cpp` for
  `src/lancet/<layer>/<unit>.h`. Multi-header units (a header that
  declares several independent classes or function families) get one
  test file per public class or function group.

- **Namespace.** Each test file's body lives inside
  `namespace lancet::<layer>::tests { ... }` — mirroring the
  `lancet::<layer>` namespace of the production code under test. The
  trailing `::tests` is mandatory; it makes "which layer does this
  test exercise?" answerable from a stack trace or symbol dump,
  prevents ODR collisions when two layers happen to name a helper the
  same (e.g. a `MakeWindow` could legitimately exist in both
  `lancet::cbdg::tests` and `lancet::core::tests`), and parallels the
  source-layer namespace one-to-one so a reader walking
  `src/lancet/<layer>/` and `tests/<layer>/` together sees the same
  `lancet::<layer>` vs. `lancet::<layer>::tests` distinction
  throughout. Inside the named namespace, unqualified lookup finds
  symbols from the parent layer — `using lancet::<layer>::Foo;` is
  redundant and should be omitted.

- File preamble: includes (project headers first via `"…"`, then
  third-party via `<…>`), then a one-line description comment if the
  file's purpose isn't obvious from its name.

- Helpers (test fixtures, factory functions, golden-data loaders)
  live at the top of the file inside the named namespace. File-private
  helpers may sit in an anonymous `namespace { … }` block *inside*
  the named namespace if internal linkage is wanted; otherwise they
  live at named-namespace scope. Cross-file shared helpers go in a
  header under `tests/<layer>/` and are declared in
  `lancet::<layer>::tests` (or `lancet::tests` if the helper genuinely
  cuts across layers — e.g. a generic temp-directory factory).

- TEST_CASEs follow, ordered to mirror the source header's foundational
  → orchestrator declaration order. The reader should be able to walk
  the test file and the source header in lockstep.

## 2. TEST_CASE naming

The TEST_CASE name describes the verified behavior in a present-tense
full sentence: subject + verb + qualifier. Examples:

- `"RevComp returns the canonical complement for each base"`
- `"OnlineStats::Merge is associative across three accumulators"`
- `"GzipOstream throws when the output path cannot be opened"`

Avoid placeholder names like `"test foo"` or `"foo works"`. The name
ends up in CI logs and in failure reports — a reader who hasn't
opened the file should be able to tell from the name alone what the
test verifies.

## 3. Tags

The required and only tag form is `[lancet][<layer>][<unit>]`. The
project does not use semantic tags such as `[property]`, `[edge]`,
`[golden]`, `[regression]`, or `[bench]`. Use Catch2's name-pattern
filtering when selective execution is needed, not extra tags.

**Per-symbol unit tag rule.** Each TEST_CASE's `<unit>` token is the
**PascalCase name of the primary public symbol that TEST_CASE
exercises**. A single test file may contain multiple unit tags when
the source header declares multiple equal-weight public symbols.

- If a TEST_CASE exercises one symbol's behavior with the others as
  setup or input, the primary symbol is the one whose behavior is
  being verified.
- If a TEST_CASE genuinely exercises a cross-cutting property over
  several equal-weight symbols at once (rare), use the
  **header-derived fallback** — the PascalCase rendering of the
  header file name.

| Source header | Public symbols | Unit tags |
|:--------------|:---------------|:----------|
| `gzip_ostream.h` | class `GzipOstream` | `[GzipOstream]` |
| `polar_coords.h` | `PolarRadius` + `PolarAngle` | `[PolarRadius]`, `[PolarAngle]` (per-symbol; `[PolarCoords]` only for cross-cutting) |
| `compute_stats.h` | class `OnlineStats` + free `Mean`/`Median`/`Minimum` | `[OnlineStats]`, `[Mean]`, `[Median]`, `[Minimum]` |
| `repeat.h` | `HammingDist` + `HasRepeat` + `HasExactRepeat` | `[HammingDist]`, `[HasRepeat]`, `[HasExactRepeat]` (per-symbol; `[Repeat]` only for cross-cutting) |

Per-symbol tags give finest-grained filterability
(`pixi run test --filter [PolarAngle]` runs only the angle tests) and
clearer failure signal (the failing tag identifies the broken
symbol). The hook in `validate_cpp_identifiers.py` polices the
`[lancet][<layer>]` portion; the unit-tag rule is enforced by this
guide and reviewed at PR time.

## 4. Test types and the quality bar

The default coverage bar (**Bar 1**) for every public function:

- One happy-path TEST_CASE.
- One edge-case TEST_CASE that hits at least one entry from the
  edge-case taxonomy (§6) for each input type.
- One error-path TEST_CASE if the function documents a failure mode
  (throws, returns `std::optional` empty, returns an error code).

Numerical or algorithm-heavy units add (**Bar 3**, where the
investment pays):

- Property test via `GENERATE` / `GENERATE_RANDOM` (or a manually
  seeded `std::mt19937_64` loop with a const literal seed) that
  documents an invariant: idempotence, monotonicity, round-trip,
  symmetry, scale invariance, permutation invariance, etc.
- Cross-validation against a reference implementation. Common
  references:
  - **An in-tree reference implementation** (e.g. a vendored C
    library compiled into the test binary alongside the C++ port).
  - **A std-lib equivalent** that the in-tree code is a tuned
    replacement for (e.g. checking a fast approximation against the
    canonical std-lib function over a parameter grid).
  - **A scientific-reference Python script** producing a precomputed
    TSV/JSON fixture (see §11 for the fixture-script layout
    convention).

Stateful units add when relevant:

- State-machine tests covering idempotent operations (Close called
  twice; Reset twice), operations after a terminal state (Write
  after Close should throw or no-op deterministically), and ordering
  violations. The standard pattern: one TEST_CASE per state
  transition, named after the transition.

## 5. SECTION vs separate TEST_CASE

- Use `SECTION` when variants share setup state. Catch2 reruns the
  TEST_CASE body from the top for each section, so common setup
  appears once and each section exercises one variant. Good for
  table-driven tests with shared scaffolding.
- Use separate TEST_CASEs when scenarios are conceptually
  independent. A SECTION-heavy TEST_CASE that mixes happy + edge +
  error paths is worse than three separate TEST_CASEs because the
  per-section reset obscures the test's intent.

## 6. Edge-case taxonomy

Pick the applicable classes per input type:

| Input type | Classes (pick at least one) |
|:-----------|:----------------------------|
| String / sequence | empty; single-char; max-realistic-length; all-same-char; mixed-case |
| Genomic sequence | also: all-N; ambiguous IUPAC codes; lowercase |
| Float | 0.0; ±0.0; ±inf; NaN; MIN_DENORM; MAX; near-boundary |
| Integer | 0; MIN; MAX; overflow boundary |
| Span / container | empty; single; SIMD-aligned; SIMD-unaligned tail |
| Group inputs (statistical) | one group empty; both empty; identical groups; asymmetric n; single-element group |
| File / IO | empty file; oversize file (compresses well, compresses poorly); post-Close write; permission-denied path |
| State machine | terminal-op called twice (idempotency); operation after terminal state; reset-then-reuse |

## 7. Random-input rules

- **Always** pin the seed with a const literal. The standard form is
  a manually seeded `std::mt19937_64` loop or
  `GENERATE(take(N, random(min, max)))` with `--rng-seed=X` set in
  the run target. **Never** use `std::random_device()` anywhere in
  `tests/`.
- The TEST_CASE name documents the property under test
  (e.g. `"<Function> is involutive on random input"`,
  `"<Algorithm> agrees with a reference implementation on random
  input"`).
- Don't use random generation for happy-path fixtures — use
  `SECTION` or hard-coded inputs. Random generation is for
  properties, not enumeration.
- Iteration count: 100–1000 per property is the default range.
  Higher if the property is cheap to evaluate, lower if it's
  expensive. Document the chosen count in a comment when it differs
  from the default.

When clang-tidy's
`bugprone-random-generator-seed` / `cert-msc32-c` / `cert-msc51-cpp`
checks flag a const-literal seed, suppress the warning at that line
with a NOLINTNEXTLINE comment that references this style guide
(`§7`). Predictability is the explicit goal in tests; the check is
designed for production code.

## 8. Assertion macros

- `REQUIRE(...)` — invariant; subsequent code in the test depends on
  it. Test aborts on failure.
- `CHECK(...)` — non-fatal; multiple checks per TEST_CASE OK; useful
  for verifying multiple independent properties without aborting.
- `REQUIRE_THROWS_AS(expr, ExceptionT)` — error-path; verifies
  exception type.
- `REQUIRE_THROWS_WITH(expr, msg)` — error-path; verifies message
  substring. Use sparingly — couples the test to message wording.
- `INFO(...)` / `CAPTURE(...)` — context for tricky assertions.
  `CAPTURE(x, y)` adds `x` and `y` to the failure log; preferred
  over hand-rolling a printf.

## 9. Determinism discipline

- All RNG seeded with a const literal at the top of the TEST_CASE
  (or in a helper function that the TEST_CASE calls with the seed
  as an argument).
- No `std::random_device()` anywhere in `tests/`.
- **No wall-clock dependencies.** Where the codebase provides a
  clock-injection seam (a function-pointer or callable parameter
  the unit-under-test accepts in its constructor), use it to inject
  a deterministic clock in the test. Tests that read `absl::Now()`
  or any other system clock are forbidden.
- **No own-filesystem dependencies outside `tests/data/`** (read-only)
  and the Catch2-managed temporary directory (writable). Where the
  unit-under-test exposes an IO sink seam (e.g. a `std::ostream&`
  constructor overload), use it to avoid temp files entirely. When
  a temp file is genuinely necessary, allocate it under
  `std::filesystem::temp_directory_path()` with a per-test-name
  scratch subdirectory and clean up at the end of the TEST_CASE.
- Tests are independent of execution order. No global mutable state
  between TEST_CASEs. Catch2 reorders tests by default; flaky tests
  whose outcome depends on prior-test side effects are bugs.

## 10. Skipping

- Use `SKIP("reason: ...");` with an explicit reason. Never `return;`
  silently from a TEST_CASE.
- Skip-on-missing-fixture is fine if the fixture is intentionally
  out-of-band (e.g. cross-validation against a multi-GB external
  dataset distributed via `data/download_test_data.sh`). Document
  *why* the fixture might be missing in the SKIP message and point
  at the regeneration command.

## 11. Test fixtures and generator scripts

Test data and the scripts that produce it live in two places:

- **Generator scripts** live in `tests/scripts/<name>.py`. One
  directory for all fixture generators across all layers, so the
  regeneration command is discoverable.
- **Generated data** lives in `tests/data/<layer>/<name>.tsv`
  (or `.json`, `.bed`, etc.). Layer-scoped because each layer's
  fixtures are consumed only by that layer's tests.

Each generator script:

- Includes a docstring explaining the regeneration command in the
  exact form a maintainer should run, e.g.:
  ```
  pixi run -e test-fixtures \
      python tests/scripts/<name>.py
  ```
- Documents the project-side semantic the fixture cross-validates
  against (the formula, the sign convention, the empty-input
  contract — whatever the C++ side asserts).
- Pins the random seeds it uses for reproducibility, so a maintainer
  regenerating on a different machine produces a byte-identical TSV.
- Pins the scientific-Python dependency versions in `pixi.toml`'s
  `test-fixtures` feature, NOT the default environment.

The C++ test reads the generated file via the
`TESTS_BASE_DATA_DIR` (or analogous per-layer) constant exposed
through `tests/test_config.h.inc`. CMake plumbs the absolute path
into the binary at configure time.

The fixture itself is committed alongside the script. Tests do not
depend on the script being runnable at test time.

## 12. Cross-references

- Layer-direction enforcement and identifier naming: `cpp_style.md`.
- Code comment conventions (the `CRITICAL:` prefix, ASCII diagrams,
  no-filler-words rule): `code_comments.md`.
- Test-data location decision tree (which fixture for which
  workflow): the `test-data-locations` skill.
- Catch2 idiom cheat sheet (deeper reference for advanced features
  — generators, sections, matchers, reporters, sharding, seed
  reproduction): the `add-cpp-test` skill's `references/`.
