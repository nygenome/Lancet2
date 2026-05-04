# C++ Code Style and Linting

The audience for this document is **anyone editing C++ in `src/`, `tests/`, or `benchmarks/`**. clang-format and clang-tidy enforce most of what follows automatically as CI gates; several rules — quote-vs-angle-bracket includes, struct member layout, the prefer-`<algorithm>`-and-Abseil rule — are enforced by code review.

Both `.clang-format` and `.clang-tidy` are extensively commented with rationale for every setting. Read them before modifying.

| Tool | Config | Enforced by |
|:-----|:-------|:------------|
| clang-format | `.clang-format` | `pixi run fmt-check` (CI gate) |
| clang-tidy | `.clang-tidy` | `pixi run lint-check` (CI gate) |

The CI workflow `lint_cpp.yml` runs both on every push and pull request. Concurrency guards (`lint-${{ github.ref }}`) cancel in-progress runs when new commits are pushed.

## Formatting rules (clang-format)

### Line width and indentation

- **100-column limit.** Wide enough for descriptive names and templates, narrow enough for side-by-side diff review.
- **2-space indent** for blocks, **4-space continuation** for wrapped args/expressions.
- **Spaces only** — no tabs.

### East const

Lancet2 enforces **east const** placement: the `const` qualifier appears *after* the type.

```cpp
// ✅ Correct (east const)
int const x = 42;
auto const& ref = container;
std::string_view const seq = read.Sequence();

// ❌ Wrong (west const)
const int x = 42;
const auto& ref = container;
```

The full qualifier order enforced by clang-format is: `inline` → `static` → `constexpr` → *type* → `const` → `volatile` → `restrict`.

### Pointer and reference alignment

Left-aligned (Google style): `Type* ptr`, `Type& ref` — not `Type *ptr` or `Type * ptr`.

### Include ordering

Includes are automatically sorted and grouped into four tiers separated by blank lines:

```cpp
// Tier 0: Main header (auto-detected from filename)
#include "lancet/cbdg/graph.h"

// Tier 1: Project headers
#include "lancet/base/types.h"
#include "lancet/cbdg/kmer.h"

// Tier 2: Third-party headers
#include "absl/container/flat_hash_map.h"
#include "spdlog/fmt/bundled/core.h"

// Tier 3: C++ standard library
#include <algorithm>
#include <vector>

// Tier 4: C standard library wrappers
#include <cmath>
#include <cstdint>
```

Do **not** manually organize includes — let clang-format handle the grouping. If `// clang-format off` is needed to protect include order (e.g., `extern "C"` blocks), document why with a one-line comment.

### Quote vs angle-bracket includes

Angle-bracket includes (`<>`) are reserved **exclusively** for C and C++ standard library headers. All other headers — project headers, third-party libraries (Abseil, spdlog, spoa, htslib, etc.) — use quote includes (`""`).

```cpp
// ✅ Correct
#include "absl/container/flat_hash_map.h"  // third-party → quotes
#include "spdlog/fmt/bundled/core.h"       // third-party → quotes
#include <algorithm>                       // C++ stdlib → angle brackets
#include <cmath>                           // C stdlib → angle brackets

// ❌ Wrong
#include <absl/container/flat_hash_map.h>  // third-party must NOT use angle brackets
```

This convention is not enforced by clang-format or clang-tidy — it must be followed manually during code review.

### Brace style

K&R "Attach" — opening brace on the same line as the statement. Empty bodies compacted to `{}`.

### Integer literal separators

Large literals get `'` separators for readability. Only applied to literals with at least 5 digits:

- Decimal, groups of 3: `1'048'576`
- Hex, groups of 2: `0xFF'00'FF`
- Binary, groups of 4: `0b1111'0000`

### Comment reflowing

clang-format **reflows long comments** to fit within the 100-column limit. This interacts directly with code comment conventions — keep comment lines concise so reflowing doesn't produce awkward breaks. Use `// clang-format off` / `// clang-format on` around ASCII art diagrams and visual tables to prevent reflowing. See `code_comments.md` for the full set of comment patterns.

### Logging and formatting

Use **spdlog macros** for all logging — `SPDLOG_INFO`, `SPDLOG_DEBUG`, `SPDLOG_WARN`, `SPDLOG_ERROR`, `SPDLOG_TRACE`, `SPDLOG_CRITICAL`. The macros are zero-cost when the corresponding log level is disabled (compile-time elided), where the function-call form is not.

Use **fmtlib** (`fmt::format`, `fmt::print`, `fmt::format_to`) for all formatting. spdlog uses fmtlib internally, so the format-string syntax is consistent across the codebase.

```cpp
// ✅ Correct
SPDLOG_INFO("Loaded {} reads in {:.2f}s", num_reads, elapsed_seconds);
auto msg = fmt::format("region {}:{}-{}", chrom, start, end);
fmt::print(stderr, "warning: {}\n", text);
```

`std::format` and `std::print` (C++20/C++23 stdlib) are **forbidden**. The at-write hook (`validate_cpp_identifiers.py`) hard-blocks both. Three reasons:

1. **Toolchain coverage.** Lancet2 supports gcc and clang at the project's pinned versions; `std::format` support has been inconsistent across compiler versions, with subtle semantic differences in formatter dispatch.
2. **Library coupling.** spdlog uses fmtlib; mixing `std::format` with spdlog logging splits the format-string ecosystem with no benefit.
3. **Locale behavior.** `std::format` is locale-sensitive by default in some implementations; fmtlib is locale-free unless explicitly opted in. For a bioinformatics tool, locale-free is the safer default.

If a third-party header pulls in `<format>` or `<print>` transitively, that is fine — the rule is about *Lancet2 source* invoking these facilities directly.

### Assertions and namespaces

**`LANCET_ASSERT` is the only assertion macro permitted in Lancet2 source.** It is defined in `src/lancet/base/assert.h` and behaves like `assert()` under `LANCET_DEBUG_MODE` builds (Debug, sanitizer, and test builds — the test target unconditionally defines `LANCET_DEBUG_MODE` so assertions fire under Release tests as well). Under Release-without-tests, `LANCET_ASSERT` is a no-op.

```cpp
#include "lancet/base/assert.h"

LANCET_ASSERT(idx < container.size());
LANCET_ASSERT(state == State::INITIALIZED && "must call Init() first");
```

**Bare `assert()` from `<cassert>` is forbidden.** The at-write hook hard-blocks it. The reason: bare `assert()` is controlled by `NDEBUG`, which Lancet2's build system flips off under Release; a bare `assert()` written in Release-deployed code is dead. `LANCET_ASSERT` makes the project's debug-mode policy explicit and consistent.

**`using namespace std` is forbidden in headers.** It pollutes every translation unit that transitively includes the header with the entire `std` namespace, breaking ADL in unpredictable ways and causing compile errors when downstream code introduces a name that collides with a `std` symbol. The at-write hook hard-blocks `using namespace std` (and `using namespace absl`) in `.h` files. In `.cpp` files, it is also strongly discouraged but not hook-blocked — prefer named `using std::vector` declarations or fully-qualified names.

```cpp
// ✅ Correct (in a .cpp file)
#include <vector>
using std::vector;
vector<int> values;

// ✅ Even better (no using-declaration)
#include <vector>
std::vector<int> values;

// ❌ Forbidden in headers; discouraged in .cpp
using namespace std;
```

## Prefer `<algorithm>` and Abseil over manual loops

Use `<algorithm>`, `<ranges>`, and Abseil utility headers (`absl/algorithm/`, `absl/strings/`, `absl/container/`) instead of hand-written loops and conditional branching wherever possible. The primary goal is **readability, maintainability, and expressiveness** — a `std::ranges::transform` or `absl::StrJoin` communicates intent more clearly than a raw `for` loop with manual index tracking.

```cpp
// ✅ Preferred — intent is immediately clear
std::ranges::transform(paths, std::back_inserter(seqs),
                       [](auto const& path) -> std::string { return std::string(path.Sequence()); });

auto const has_case = std::ranges::any_of(samps, [](auto const& s) {
  return s.TagKind() == cbdg::Label::CASE;
});

auto const total = std::accumulate(components.begin(), components.end(), u64{0},
                                   [](u64 sum, auto const& comp) -> u64 { return sum + comp.size() - 1; });

// ❌ Avoid — manual loop obscures intent
std::vector<std::string> seqs;
for (auto const& path : paths) {
  seqs.push_back(std::string(path.Sequence()));
}
```

**Exception:** Do not use `<algorithm>` or Abseil when it would hurt readability, unnecessarily complicate the code, or degrade performance in any meaningful way. A simple 3-line loop that is immediately clear can be better than a contorted `std::transform` with multiple lambdas. Use judgment — the goal is clarity, not dogma.

## Static analysis (clang-tidy)

### Naming conventions

Enforced by `readability-identifier-naming` checks:

| Element | Convention | Example |
|:--------|:-----------|:--------|
| Classes / Structs | `PascalCase` | `VariantCall`, `AsyncWorker` |
| Methods / Functions | `PascalCase` | `BuildSequence`, `FindEdges` |
| Member Variables | `mPascalCase` | `mKmer`, `mEdges`, `mCurrK` |
| Local Variables | `snake_case` | `new_edge`, `comp_id` |
| Parameters | `snake_case` | `stop_token`, `window_length` |
| Constants (`constexpr` / `static const`) | `UPPER_CASE` | `DEFAULT_MIN_KMER_LEN` |
| Enum Values | `UPPER_CASE` | `PLUS_PLUS`, `FWD`, `UNKNOWN` |
| Type Aliases | `PascalCase` | `NodeID`, `EdgeList` |
| Template Parameters | `PascalCase` | `Args`, `TagResultValue` |
| Namespaces | `snake_case` | `lancet`, `cbdg`, `hts` |
| Macros | `UPPER_CASE` | `LANCET_DEVELOP_MODE` |

**Ignored patterns:** type aliases `i8`, `u32`, `f64`, `usize` are exempted. Include guard macros (`SRC_LANCET_*_H_`) are exempted. STL iterator methods (`begin`, `end`, `cbegin`, `cend`) are exempted.

### Struct/class memory layout

All struct and class member variables **must** be declared in descending order of alignment size to minimize compiler-inserted padding. This is not enforced by clang-format or clang-tidy — it is enforced in code review.

**Alignment tiers** (largest → smallest):

| Tier | Size | Types |
|:-----|:-----|:------|
| 8B | 8 bytes | `i64`, `u64`, `usize`, `f64`, pointers (`T*`), `std::string`, `std::vector<T>`, `absl::flat_hash_map`, `std::unique_ptr<T>`, `absl::InlinedVector<T,N>` |
| 4B | 4 bytes | `i32`, `u32`, `f32`, `std::array<T,N>` where T is 4B |
| 2B | 2 bytes | `i16`, `u16` |
| 1B | 1 byte | `i8`, `u8`, `bool`, `enum class : i8` |

**Constructor initializer lists** must match declaration order. The compiler warns (`-Wreorder`) if the initializer list order doesn't match — always write initializer lists in declaration order.

```cpp
// ✅ Correct — declaration and init list both follow 8B → 4B → 1B
struct ReadEvidence {
  f64 mScore;       // 8B
  i32 mCount;       // 4B
  u8 mQuality;      // 1B
  bool mIsForward;  // 1B

  ReadEvidence(f64 s, i32 c, u8 q, bool f)
      : mScore(s), mCount(c), mQuality(q), mIsForward(f) {}
};

// ❌ Wrong — bool before f64 wastes 7 bytes of padding
struct BadLayout {
  bool mIsForward;  // 1B + 7B padding before f64
  f64 mScore;       // 8B
  u8 mQuality;      // 1B + 3B padding before i32
  i32 mCount;       // 4B
};
```

When using designated initializers or named field assignment (e.g., `result.mScore = x`), order doesn't matter for correctness, but **declaration order still determines padding** — so the declaration must follow the tier ordering.

The size annotations (`// 8B`, `// 4B`, etc.) on member variable declarations serve double duty: they document the field and make alignment compliance visible at a glance. The comment convention is documented in `code_comments.md` § "Member Variable Comments"; the layout rule itself is what this section enforces.

### Check categories

The following check families are enabled, with curated exclusions documented in `.clang-tidy`:

| Category | Coverage | Key Exclusions |
|:---------|:---------|:---------------|
| `bugprone-*` | Bug detection | `-easily-swappable-parameters` (too noisy for graph APIs) |
| `cert-*` | CERT secure coding | `-err58-cpp` (static constexpr never throws) |
| `cppcoreguidelines-*` | Core Guidelines | FFI checks disabled (htslib/minimap2 require raw arrays, varargs) |
| `modernize-*` | C++20 modernization | `-use-trailing-return-type` (breaks lambda deduction in GCC), `-use-ranges` (breaks types without `<=>`) |
| `performance-*` | Performance | All enabled; `string_view` and `Span` exempted from value-param check |
| `readability-*` | Readability and naming | `-magic-numbers` (scientific constants), `-redundant-member-init` (explicit init is convention) |
| `misc-*` | Miscellaneous | `-no-recursion` (graph DFS/BFS is inherently recursive), `-include-cleaner` (breaks include order) |

### Warnings promoted to errors

These checks are hard CI errors — they will fail the build:

- All `bugprone-*` (except `exception-escape`)
- All `cert-*`, `performance-*`, `portability-*`
- `modernize-use-override`, `modernize-use-nullptr`, `modernize-use-emplace`
- `readability-container-contains`, `readability-container-size-empty`
- `misc-const-correctness`

### Function complexity thresholds

| Metric | Threshold |
|:-------|:----------|
| Cognitive complexity | 35 |
| Statement count | 150 |
| Line count | 200 |
| Parameter count | 6 |
| Nesting depth | 4 |

Functions exceeding these thresholds should be decomposed. If a function legitimately needs to exceed a threshold (e.g., a CIGAR walk), use a scoped `// NOLINTNEXTLINE(readability-function-cognitive-complexity)` with a rationale comment.

### NOLINT usage rules

Suppression comments are acceptable **only** when:

1. The code is correct and the check is a false positive (e.g., C FFI interop).
2. The suppression uses scoped form: `NOLINTNEXTLINE(check-name)` or a `NOLINTBEGIN(check-name)` / `NOLINTEND(check-name)` block — **never** inline `NOLINT` (it pollutes the statement line) and **never** bare/unscoped (omitting the check name).
3. A rationale comment explains *why* the suppression is needed.

```cpp
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast) -- htslib C FFI requires raw cast
auto const* aptr = reinterpret_cast<unsigned long long const*>(first.data());
```

`NOLINTBEGIN` / `NOLINTEND` blocks are acceptable for groups of related suppressions (e.g., struct members that must be `const&`):

```cpp
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
absl::Span<u8 const> const mQuery;         // 16B
absl::Span<u8 const> const mTarget;        // 16B
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
```

The first answer to "should I suppress this?" is always "no — fix the underlying issue." Suppressions are for cases where the check is genuinely wrong about this code, not for cases where fixing the check finding is inconvenient.

## Running locally

All build tools, linters, and documentation generators are managed by [pixi](https://pixi.sh) and installed into the project-local `.pixi/` directory. Every developer operation is available as a `pixi run <task>` command. Run `pixi task list` for the full catalog.

### Build and test

```bash
# Configure and build in Release mode
pixi run configure-release            # CMake Release (static, no cloud I/O)
pixi run build-release                # Build Release (auto-runs configure-release first)

# Configure and build in Debug mode (with tests)
pixi run configure-debug              # CMake Debug (tests ON, compile_commands.json)
pixi run build-debug                  # Build Debug (auto-runs configure-debug first)

# Run tests (auto-builds Release first; LANCET_DEBUG_MODE is defined on the
# test target unconditionally so LANCET_ASSERT keeps firing under Release)
pixi run test

# Run benchmarks in Release mode
# Step 1: configure with benchmarks enabled (one-time)
pixi run configure-release -- -DLANCET_BENCHMARKS=ON
# Step 2: build and run
pixi run bench                        # builds Release, then runs benchmarks
```

Both `configure-release` and `configure-debug` accept extra CMake flags via `--`. The preset flags (`BUILD_TYPE`, `BUILD_STATIC`, `ENABLE_CLOUD_IO`) define the build "personality." Optional flags (`BENCHMARKS`, `PROFILE_MODE`, `BUILD_ARCH`) persist in CMake's cache once set:

```bash
# Release + cloud I/O + custom CPU arch
pixi run configure-release -- -DLANCET_ENABLE_CLOUD_IO=ON -DLANCET_BUILD_ARCH=native

# Debug + benchmarks (e.g. for linting all targets)
pixi run configure-debug -- -DLANCET_BENCHMARKS=ON
```

### Formatting

```bash
pixi run fmt-check                    # dry-run — reports files that need formatting
pixi run fmt-fix                      # applies formatting in-place to all sources
```

For file-specific or directory-specific formatting, invoke the script directly:

```bash
python3 scripts/run_clang_format.py src/lancet/cbdg
python3 scripts/run_clang_format.py --fix src/lancet/caller/genotyper.cpp
```

### Static analysis

```bash
pixi run lint-check                   # check-only (auto-builds Release first for compile_commands.json)
pixi run lint-all                     # run fmt-check + lint-check + iwyu-check
```

> [!CAUTION]
> **Clang-tidy auto-fix is not supported in this project.** The `lint-fix` pixi task and the `--fix` flag on `scripts/run_clang_tidy.py` were removed because clang-tidy auto-fixes have historically broken compilation and produced unreadable code (modernize-use-trailing-return-type rewriting lambdas, modernize-use-ranges breaking sort, misc-include-cleaner reordering includes). Always resolve clang-tidy warnings manually by reading the diagnostic, understanding the root cause, and writing the fix yourself.

For custom build directories, invoke the script directly:

```bash
python3 scripts/run_clang_tidy.py --build-dir cmake-build-release
```

### Documentation

```bash
pixi run docs-serve                   # local preview at http://127.0.0.1:8000
pixi run docs-build                   # strict validation (zero-warning required)
pixi run docs-export                  # export all docs to single Markdown file
```

### Release and deployment

```bash
pixi run version-bump                 # bump patch version, update changelog, push
pixi run version-bump --kind minor    # bump minor version
pixi run version-bump --kind major    # bump major version
pixi run update-changelog             # regenerate changelog standalone (requires Go)
pixi run conda-build-local            # build local conda package via rattler-build
pixi run docker-push                  # build and push Docker image to GCR
```

### CI behaviour

The `lint_cpp.yml` workflow runs both linting checks on every push and pull request — the same tools invoked by `pixi run fmt-check` and `pixi run lint-check`:

1. **clang-format** — dry-run check on all `src/lancet/`, `tests/`, `benchmarks/` files.
2. **clang-tidy** — full build + analysis with `compile_commands.json`.

Concurrency guards (`lint-${{ github.ref }}`) cancel in-progress runs when new commits are pushed.

## Git commit messages

Lancet2 uses [Conventional Commits](https://www.conventionalcommits.org/) with a simplified pattern — no scopes. The changelog (`CHANGELOG.md`) is auto-generated from commit messages by `git-chglog` (via `pixi run version-bump` or `pixi run update-changelog`). Only commits matching recognized types appear in the changelog; all other commits are silently excluded.

### Format

```
type: concise imperative summary of the change

Optional body — explains what changed and why. Wrap all lines at 72
characters. Use the same language standards as website documentation:
active voice, no filler, no hedging, quantify when possible.
```

The commit message is documentation for someone trying to understand what changed *without reading the diff*. Write it for that reader.

- **Subject line** (first line): standalone summary of the entire commit. Must be self-contained — a developer scanning `git log --oneline` should understand the change from this line alone. Lowercase, imperative mood ("add" not "added"), no trailing period.
- **Blank line**: separates subject from body. Required if a body is present.
- **Body** (optional): explains *what* changed and *why*. Wrap every line at **72 characters** — this ensures readability in `git log`, terminal pagers, and email-based review tools. Use bullet points or short paragraphs for multi-part changes. Do not restate the diff — describe the intent, rationale, and anything non-obvious.

**Type** is one of exactly four values:

| Type | Changelog Section | When to Use |
|:-----|:------------------|:------------|
| `feat` | New Features | New user-facing functionality: CLI flags, VCF fields, algorithm capabilities |
| `fix` | Bug Fixes | Corrects incorrect behaviour: wrong output, crashes, logic errors |
| `perf` | Performance Improvements | Measurable performance change: reduced allocations, faster hot paths, lower memory |
| `chore` | Refactoring | Everything else: refactoring, dependency updates, CI changes, docs, formatting, tooling |

### Examples

```bash
# ✅ Subject-only (most commits) — self-contained one-liner
feat: add AHDD and PDCV artifact metrics to VCF output
fix: correct NM spec violation in ASMD computation
perf: eliminate per-read heap allocations via zero-copy bam1_t proxying
chore: resolve clang-tidy warnings in tests
```

```
# ✅ Subject + body — for changes that benefit from context
feat: add Dirichlet-Multinomial genotype likelihoods and CMLOD field

Replace the per-allele binomial model with a Dirichlet-Multinomial
that accounts for overdispersion from PCR duplicates. The collapsed
model log-odds (CMLOD) is emitted as a new Float FORMAT field.

- VariantBuilder: thread alpha prior through GenotypeLikelihoods()
- variant_call.h: add mCMLOD member (8B, f64)
- VCF header: register CMLOD with Number=1, Type=Float
```

```bash
# ❌ Bad — wrong type, past tense, vague, or missing type
Added new metrics                      # missing type, past tense
feat: Update stuff                     # uppercase subject, vague
refactor: rename function              # "refactor" is not a recognized type, use "chore"
fix: fixes the bug                     # not imperative ("fix" not "fixes")
```

### Breaking changes

For changes that break backward compatibility (e.g., renamed CLI flags, removed VCF fields, changed output format), add a `BREAKING CHANGE` paragraph to the commit body:

```
feat: replace STR_SCORE with LCR_SCORE in VCF output

BREAKING CHANGE: The STR_SCORE FORMAT field has been removed and replaced
by LCR_SCORE. Downstream parsers that depend on STR_SCORE will need to be
updated to read LCR_SCORE instead.
```

Breaking changes produce a dedicated "BREAKING CHANGE" section in the changelog.

### What NOT to do

- **Do not use scopes.** The `.chglog` parser does not extract scopes (`feat(graph): ...` is parsed but the scope is discarded).
- **Do not capitalize the subject.** The changelog renders subjects verbatim; capitalized subjects look inconsistent.
- **Do not use types other than `feat` / `fix` / `perf` / `chore`.** Anything else is silently excluded from the changelog.
- **Do not end the subject with a period.** It reads as a sentence fragment in the changelog bullet list.

## Cross-references

- The comment conventions that interact with comment-reflowing rules and the size-annotation convention live in `code_comments.md`.
- The conventions for adding or editing tests under `tests/` (per-symbol PascalCase tag rule, `namespace lancet::<layer>::tests`, determinism discipline, property-test rules, test-fixture-script layout) live in `test_style.md`.
- The synchronization rule for keeping code consistent with website docs and dev docs lives in `sync_and_verification.md`.
- The voice and tone rules for the user-facing website docs are in `website_docs.md` — different audience, different bar.
