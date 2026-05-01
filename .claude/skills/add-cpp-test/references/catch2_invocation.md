# Running the test binary flexibly

This reference is loaded on demand by the `add-cpp-test` skill when
running tests directly (instead of via `pixi run test`). The binary
is `cmake-build-release/tests/TestLancet2`; all Catch2 flags apply.


`pixi run test` builds Release and runs `cmake-build-release/tests/TestLancet2` with no args, which runs all non-hidden tests in random order. For day-to-day work that's the right invocation. For iterative debugging, drive the binary directly with the Catch2 CLI — the binary is the same executable, so all Catch2 flags apply.

### Selecting which tests to run

The first positional arg is a test spec. Without one, all non-hidden tests run.

```bash
# Build first if needed.
pixi run --quiet build-release

# Run by exact name.
./cmake-build-release/tests/TestLancet2 "graph builds correctly"

# Run by wildcard (escape * in the shell).
./cmake-build-release/tests/TestLancet2 "graph*"

# Run all tests in a layer (tag selection).
./cmake-build-release/tests/TestLancet2 "[cbdg]"

# Multiple tags AND-ed: tests that have both tags.
./cmake-build-release/tests/TestLancet2 "[cbdg][hash]"

# Multiple tags OR-ed: tests that have either tag (comma).
./cmake-build-release/tests/TestLancet2 "[cbdg],[hts]"

# Negate a tag: skip slow tests.
./cmake-build-release/tests/TestLancet2 "~[.slow]"

# Combine: caller layer except slow.
./cmake-build-release/tests/TestLancet2 "[caller]~[.slow]"

# Run a hidden test by addressing it explicitly.
./cmake-build-release/tests/TestLancet2 "[.slow]"
```

### Running a specific section

When a test case has multiple sections and you want one, use `-c` / `--section`:

```bash
./cmake-build-release/tests/TestLancet2 "vector reset behaviour" -c "clear() empties the vector"
```

For nested sections, repeat `-c` for each level. The flag stacks.

### Running a specific generator value

When debugging a `GENERATE`-driven test, isolate one value with `-g` / `--generator-index`:

```bash
./cmake-build-release/tests/TestLancet2 "Phred conversion round-trips" -g 2
```

The index is zero-based.

### Discovery: list tests, tags, and reporters

```bash
./cmake-build-release/tests/TestLancet2 --list-tests           # all test names + tags
./cmake-build-release/tests/TestLancet2 --list-tests "[cbdg]"  # filtered
./cmake-build-release/tests/TestLancet2 --list-tags            # all tags + counts
./cmake-build-release/tests/TestLancet2 --list-reporters       # available output formats
./cmake-build-release/tests/TestLancet2 --verbosity quiet --list-tests   # names only

# Pipe filtered names into a file for batch re-runs.
./cmake-build-release/tests/TestLancet2 --list-tests --verbosity quiet "[caller]" > caller_tests.txt
./cmake-build-release/tests/TestLancet2 -f caller_tests.txt
```

### Reproducing a failed random-order run

The test binary defaults to random order. When a test fails intermittently and order matters, the `Randomness seeded to: <N>` line in the output identifies the seed:

```bash
./cmake-build-release/tests/TestLancet2 --rng-seed 0xCAFEBABE
./cmake-build-release/tests/TestLancet2 --order decl    # declaration order, deterministic
./cmake-build-release/tests/TestLancet2 --order lex     # lexicographic by name
```

Random order is subset-stable: filtering doesn't change relative test order at a fixed seed. Use this to bisect interference: if `--rng-seed N "[cbdg]"` reproduces and `--rng-seed N "[cbdg][hash]"` does not, the failing test depends on something earlier in the layer.

### Debugging modes

```bash
# Show output for passing tests too (default suppresses).
./cmake-build-release/tests/TestLancet2 -s "graph builds correctly"

# Break into the debugger on first failure.
./cmake-build-release/tests/TestLancet2 -b "graph builds correctly"

# Stop after first failure (any kind).
./cmake-build-release/tests/TestLancet2 -a

# Stop after N assertion failures.
./cmake-build-release/tests/TestLancet2 -x 5

# Make whitespace visible in string-comparison failures.
./cmake-build-release/tests/TestLancet2 -i

# Report per-test wall time.
./cmake-build-release/tests/TestLancet2 -d yes
./cmake-build-release/tests/TestLancet2 -D 0.5    # only tests over 0.5s

# Treat as warnings: tests with no assertions, unmatched specs, infinite generators.
./cmake-build-release/tests/TestLancet2 -w NoAssertions -w UnmatchedTestSpec
```

The `NoAssertions` warning is worth using on every CI run — it catches accidentally-empty test bodies and unreached leaf sections.

### Reporters

The default Console reporter is human-friendly. For machine-readable output (CI, regression tracking), pick another reporter via `-r`:

```bash
./cmake-build-release/tests/TestLancet2 -r JUnit -o results.xml
./cmake-build-release/tests/TestLancet2 -r XML -o results.xml
./cmake-build-release/tests/TestLancet2 -r SonarQube -o results.xml
./cmake-build-release/tests/TestLancet2 -r compact          # one-line per failure
```

Multiple reporters in a single run write to different sinks:

```bash
./cmake-build-release/tests/TestLancet2 \
    -r JUnit::out=results.xml \
    -r console::out=-::colour-mode=ansi
```

That writes JUnit XML to disk and human-friendly output to stdout simultaneously. Useful in CI where you want both.

### Sharding (parallelism)

Catch2 splits test cases evenly into N shards, identified by index 0..N-1. Naïvely combining `--shard-count` with `--order rand` breaks coverage because each shard would have its own seed; share the seed across shards explicitly:

```bash
SEED=0xBEEF
./cmake-build-release/tests/TestLancet2 --order rand --rng-seed $SEED --shard-count 3 --shard-index 0 &
./cmake-build-release/tests/TestLancet2 --order rand --rng-seed $SEED --shard-count 3 --shard-index 1 &
./cmake-build-release/tests/TestLancet2 --order rand --rng-seed $SEED --shard-count 3 --shard-index 2 &
wait
```

Use shard count > core count to avoid one long-tailed test bottlenecking a shard.

### Recommended invocations

| Scenario | Invocation |
|:---|:---|
| Default day-to-day | `pixi run test` |
| Iterating on a single test | `./cmake-build-release/tests/TestLancet2 "test name"` |
| Iterating on a layer | `./cmake-build-release/tests/TestLancet2 "[caller]~[.slow]"` |
| Reproducing a random-order failure | `./cmake-build-release/tests/TestLancet2 --rng-seed <reported-seed>` |
| Pre-merge sanity (exhaustive, with no-assertion check) | `./cmake-build-release/tests/TestLancet2 --order rand --warn NoAssertions` |
| CI run with machine-readable output | `./cmake-build-release/tests/TestLancet2 -r JUnit::out=results.xml -r console::out=- --warn NoAssertions` |
| Slow tests only (normally hidden) | `./cmake-build-release/tests/TestLancet2 "[.slow]"` |

The `pixi run test` task is fine when the wrapper's `depends-on` build step is what you want; drive the binary directly when you need any of the flexibility above.
