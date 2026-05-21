# Catch2 v3.14.0 idiom cheat sheet

This reference is loaded on demand by the `add-cpp-test` skill when
writing tests. It documents the Catch2 features Lancet2 actually uses
(or could use) and surfaces the parts most often misremembered.


This section is the reference for the Catch2 features Lancet2 actually uses or could use. The Catch2 docs live at [github.com/catchorg/Catch2/tree/v3.14.0/docs](https://github.com/catchorg/Catch2/tree/v3.14.0/docs); this cheat sheet surfaces the parts most often misremembered or under-used. Read it when you find yourself reaching for repetitive scaffolding that a Catch2 feature would replace.

### Assertion families

The two-axis assertion model: `REQUIRE` aborts the test case on failure; `CHECK` records and continues. Use `REQUIRE` when subsequent assertions depend on this one being true; use `CHECK` when assertions are independent and you want to see all failures from a run.

| Form | Semantics |
|:---|:---|
| `REQUIRE(expr)` / `CHECK(expr)` | Truthy. Catch2 decomposes the expression and reports both sides. |
| `REQUIRE_FALSE(expr)` / `CHECK_FALSE(expr)` | Falsy. Use this rather than `REQUIRE(!expr)` — `!` cannot be decomposed. |
| `REQUIRE_NOTHROW(expr)` / `CHECK_NOTHROW(expr)` | Expression must not throw. |
| `REQUIRE_THROWS(expr)` / `CHECK_THROWS(expr)` | Expression must throw any exception. |
| `REQUIRE_THROWS_AS(expr, ExceptionType)` | Expression must throw an exception of exactly this type. The macro adds `const&` itself; do not write it. |
| `REQUIRE_THROWS_WITH(expr, "msg" or matcher)` | Expression must throw, and the exception's `what()` must match. |
| `REQUIRE_THROWS_MATCHES(expr, ExceptionType, matcher)` | Combined type + content match. |
| `REQUIRE_THAT(value, matcher)` / `CHECK_THAT(value, matcher)` | Matcher-based assertion. See "Matchers" below. |

Two limitations to remember. First, `&&` and `||` inside an assertion expression cannot be decomposed; either parenthesize the whole expression or split into two assertions. Second, multi-argument macros (`REQUIRE_THROWS_AS`, `TEST_CASE_METHOD` with templated fixtures) choke on commas inside template parameters: use a typedef (`using Pair = std::pair<int, int>;`) or extra parentheses.

### TEST_CASE and SECTION

`TEST_CASE("name", "[tags]")` is the basic test container. The combination of name and tags must be unique within the executable. `SECTION("name")` nests inside a `TEST_CASE`; each leaf section runs in its own pass through the test case, with the test case body re-executing from the top each time. This makes setup-then-multiple-checks workflows clean without traditional fixture overhead:

```cpp
TEST_CASE("vector reset behaviour", "[base]") {
    std::vector<int> v{1, 2, 3};
    REQUIRE(v.size() == 3);

    SECTION("clear() empties the vector") {
        v.clear();
        REQUIRE(v.empty());
    }
    SECTION("resize(0) empties the vector") {
        v.resize(0);
        REQUIRE(v.empty());
    }
}
```

Both sections run, both start fresh from the `v{1, 2, 3}` initialization. Sections nest arbitrarily; the test case re-runs once per leaf section.

### Tags

Tags appear in the second arg of `TEST_CASE` as bracketed strings: `"[caller][slow]"`. Multiple tags are concatenated. Special tags begin with non-alphanumerics:

| Tag | Behaviour |
|:---|:---|
| `[.]` or `[.foo]` | Hidden by default; runs only when explicitly selected by name or tag. Useful for slow/integration tests. |
| `[!throws]` | Test is expected to throw; excluded when `--nothrow` is set. |
| `[!mayfail]` | Failures don't fail the test (still reported). For known-broken tests. |
| `[!shouldfail]` | Inverted: passing fails the test. For tracking accidental fixes. |
| `[!benchmark]` | Hidden by default; this test contains `BENCHMARK` blocks. |
| `[!nonportable]` | Documentation only; signals the test depends on platform behaviour. |

Lancet2 convention: tag every test with its layer (`[base]`, `[hts]`, `[cbdg]`, `[caller]`, etc.). Slow tests get `[.slow]` so they're hidden from the default run. Tests that exercise specific subsystems can carry component tags (`[probe]`, `[graph]`, `[scorer]`).

### Type-parametrised tests

When the same logic should be exercised over multiple types, `TEMPLATE_TEST_CASE` runs the body once per type with `TestType` aliased to each in turn:

```cpp
#include <catch2/catch_template_test_macros.hpp>

TEMPLATE_TEST_CASE("integer hashing is collision-free over small range",
                   "[base][hash]", u32, u64) {
    std::set<TestType> seen;
    for (TestType i = 0; i < 1024; ++i) {
        seen.insert(my_hash(i));
    }
    REQUIRE(seen.size() == 1024);
}
```

Variants: `TEMPLATE_PRODUCT_TEST_CASE` (cartesian product of template-types and template-args), `TEMPLATE_LIST_TEST_CASE` (types from a `std::tuple` for reuse across test cases). For non-type template parameters, `TEMPLATE_TEST_CASE_SIG` accepts a signature.

Wrap multi-parameter types in extra parentheses: `(std::map<int, std::string>)`. Otherwise the preprocessor counts the comma as an argument separator.

### Generators

`GENERATE` produces a sequence of values; the test case runs once per generated value, like an outer SECTION:

```cpp
#include <catch2/generators/catch_generators.hpp>

TEST_CASE("Phred conversion round-trips", "[base][phred]") {
    auto q = GENERATE(0, 10, 20, 30, 40, 60);
    auto p = phred_to_prob(q);
    REQUIRE(prob_to_phred(p) == q);
}
```

Multiple `GENERATE` calls in the same scope produce a cartesian product. Useful built-in generators: `range(start, end)`, `range(start, end, step)`, `take(n, gen)`, `random(lo, hi)`, `chunk(n, gen)`, `filter(predicate, gen)`, `from_range(begin, end)`. For complex test data, `GENERATE_REF` and `GENERATE_COPY` control capture semantics. Use `--generator-index` on the command line to run one specific value during debugging.

### Fixtures

`TEST_CASE_METHOD(Fixture, "name", "[tags]")` creates a fresh fixture instance per test case (and per leaf section). Use this when the setup is expensive enough to warrant a class but cheap enough to recreate per test:

```cpp
class GraphTestFixture {
protected:
    cbdg::Graph mGraph;
public:
    GraphTestFixture() : mGraph(make_test_graph()) {}
};

TEST_CASE_METHOD(GraphTestFixture, "graph builds correctly", "[cbdg]") {
    REQUIRE(mGraph.NodeCount() > 0);
}
```

For fixtures that should persist across all test cases of a class (lazy initialization, expensive setup), use `TEST_CASE_PERSISTENT_FIXTURE`. Templated fixtures: `TEMPLATE_TEST_CASE_METHOD`. With templated types whose parameter list contains commas, wrap the whole class name in parentheses: `TEST_CASE_METHOD((Fixture<int, int>), ...)`.

### Matchers

Matchers compose with `&&`, `||`, `!`. Useful built-ins:

```cpp
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::EndsWith;
using Catch::Matchers::Equals;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

REQUIRE_THAT(some_string, ContainsSubstring("oops") && EndsWith(".txt"));
REQUIRE_THAT(some_vector, Equals(std::vector<int>{1, 2, 3}));
```

For floating-point comparisons (which Lancet2 does heavily — genotype likelihoods, scoring, posterior base quality), prefer the matchers over raw `==`. Three matchers cover the spectrum:

- `WithinAbs(target, margin)` — absolute tolerance. Use when the target's scale is known and bounded.
- `WithinRel(target, eps)` — relative tolerance. Default `eps` is `numeric_limits<FP>::epsilon * 100`. Use for ratios or when comparing values whose magnitude varies.
- `WithinULP(target, maxUlpDiff)` — units in the last place. Most precise for values that are computed from an arithmetic chain expected to be bit-exact.

Combine for robustness: a value should be within relative tolerance for non-zero targets but within absolute tolerance near zero (where relative tolerance breaks down):

```cpp
REQUIRE_THAT(scoring_output,
    WithinRel(expected, 0.001) || WithinAbs(0, 1e-9));
```

### Logging

Inside a `TEST_CASE`, the `INFO("text" << values)` macro attaches a stream-style message that prints only when a subsequent assertion fails. `CAPTURE(var1, var2)` is shorthand for printing names and values of variables. Both are scoped — they apply to assertions within the enclosing scope.

```cpp
TEST_CASE("k-mer hashing", "[cbdg]") {
    auto k = GENERATE(15, 21, 31);
    INFO("k = " << k);
    auto kmer = make_test_kmer(k);
    CAPTURE(kmer.AsString());
    REQUIRE(kmer.Hash() != 0);
}
```

`UNSCOPED_INFO` and `UNSCOPED_CAPTURE` (added in 2.7.0 / 3.13.0) attach to the next assertion only and are not scope-tied — useful inside helper functions:

```cpp
void verify_invariants(const Graph& g) {
    UNSCOPED_INFO("node count = " << g.NodeCount());
    REQUIRE(g.IsValid());
}
```

### Skipping at runtime

`SKIP("reason")` (added in 3.3.0) skips a test case at runtime — neither pass nor fail. Use when prerequisites are missing (test data file absent, hardware capability not present) rather than failing the test for unrelated reasons:

```cpp
TEST_CASE("HCC1395 truth comparison", "[caller][.slow]") {
    if (!std::filesystem::exists(LANCET_FULL_DATA_DIR "/HCC1395_truth.vcf.gz")) {
        SKIP("HCC1395 truth VCF not present in data/");
    }
    // ...
}
```

A pre-`SKIP` `CHECK`/`REQUIRE` failure still fails the test. `FAIL("reason")`, `FAIL_CHECK("reason")`, and `SUCCEED("reason")` are explicit terminators.

### Benchmarks (live in `benchmarks/`, not `tests/`)

Lancet2's benchmarks live in `benchmarks/` and use **google/benchmark**, not Catch2's `BENCHMARK` macro. The two frameworks are unrelated; do not confuse them. The `profile-and-optimize` skill covers the benchmark workflow and points at its own `references/google_benchmark_idioms.md` for the macro surface (`BENCHMARK`, `BENCHMARK_F`, `DoNotOptimize`, `ClobberMemory`, fixtures, ranges). This skill (`add-cpp-test`) is for unit/integration tests in `tests/` only.
