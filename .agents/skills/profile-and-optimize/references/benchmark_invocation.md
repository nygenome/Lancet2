# Benchmark binary invocation reference

This reference is loaded on demand by the `profile-and-optimize` skill
when running the benchmark binary directly (instead of via the default
`pixi run --quiet bench`). It documents the flag surface, common patterns
for variance reduction, and the recipe table at the end.


`pixi run --quiet bench` builds Release and runs `cmake-build-release/benchmarks/BenchmarkLancet2` with no args, which runs every registered benchmark once. For day-to-day cross-validation that's the right invocation. For iterative measurement work, drive the binary directly — every `--benchmark_*` flag has a matching `BENCHMARK_*` environment variable, with the flag winning if both are set.

### Selecting a subset

`--benchmark_filter=<regex>` selects which benchmarks to run by display name (which includes args, like `BM_HashKmer/21`):

```bash
pixi run --quiet build-release

# Filter by full name
./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_filter=BM_HashKmer

# Filter to one parameterization
./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_filter='BM_HashKmer/21'

# Anchor to start (regex semantics — escape special chars)
./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_filter='^BM_Hash'

# Multiple via alternation
./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_filter='BM_HashKmer|BM_RepeatDetect'
```

`--benchmark_list_tests=true` enumerates registered benchmarks without running them — useful for piping into a filter file.

### Stability and statistics

A single run is noisy. The library's built-in remedy is repetition with statistical aggregation:

```bash
# Run each benchmark 10 times, report mean/median/stddev/CV.
./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_repetitions=10

# Show only the aggregates (less console noise).
./cmake-build-release/benchmarks/BenchmarkLancet2 \
    --benchmark_repetitions=10 \
    --benchmark_report_aggregates_only=true

# Run each benchmark for at least 2 seconds before reporting.
./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_min_time=2.0s

# Warm up before measuring (per-iteration warm-up).
./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_min_warmup_time=1.0
```

The coefficient of variation (CV) printed alongside the mean is the right number to read. Anything above ~10% means the measurement environment is too noisy to trust the comparison; pin the CPU governor, isolate cores, or increase iterations.

### Output format and files

```bash
# Console (default), JSON, or CSV.
./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_format=json

# Write structured output to a file (separate from console).
./cmake-build-release/benchmarks/BenchmarkLancet2 \
    --benchmark_out=/tmp/bench.json \
    --benchmark_out_format=json

# Display custom counters as their own columns instead of inline.
./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_counters_tabular=true
```

Lancet2 already commits some baseline JSON files under `benchmarks/results/` (e.g., `extractor_bench.AMD_Milan_results.json`); use the matching `--benchmark_format=json --benchmark_out=...` invocation when capturing a new baseline so the file format is consistent.

### Dry runs and other utilities

```bash
# Single-iteration sanity check (compile-and-link smoke test).
./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_dry_run

# Globally set the time unit shown in the report.
./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_time_unit=us

# Add a key/value to the context block (e.g., commit hash).
./cmake-build-release/benchmarks/BenchmarkLancet2 \
    --benchmark_context=commit=$(git rev-parse --short HEAD)
```

### Comparing runs with `compare.py`

Google/benchmark ships a `compare.py` script (`tools/compare.py` in the source tree). It runs a Mann-Whitney U-test on two JSON outputs and reports per-benchmark deltas with statistical significance. The basic invocation:

```bash
./compare.py benchmarks /path/to/before.json /path/to/after.json
./compare.py benchmarks ./before-binary ./after-binary           # runs them first
./compare.py filters    ./binary 'BM_OldImpl' 'BM_NewImpl'       # same binary, two filters
```

Output columns are `Time`, `CPU`, `Time Old/New`, `CPU Old/New`. Negative deltas are speedups; the U-test p-value (when available) tells you whether the delta is statistically meaningful or noise. The script depends on `scipy`; install it in the relevant pixi environment if you reach for it.

### Reducing variance

If your CV is consistently above 10%, the problem is the environment, not the benchmark. The library docs list the practical mitigations on Linux:

```bash
# Switch to performance governor (most impactful single change).
sudo cpupower frequency-set --governor performance

# Disable Turbo Boost (reduces frequency variance).
echo 0 | sudo tee /sys/devices/system/cpu/cpufreq/boost

# Pin to one core.
taskset -c 0 ./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_filter=...

# Raise scheduling priority.
sudo nice -n -20 ./cmake-build-release/benchmarks/BenchmarkLancet2 --benchmark_filter=...
```

The library also emits a `***WARNING*** CPU scaling is enabled` line when it detects the issue. If you see it, fix the environment before trusting any number.

### Recommended invocations

| Scenario | Invocation |
|:---|:---|
| Default cross-validation | `pixi run --quiet bench` |
| One benchmark, repeated for stability | `--benchmark_filter=NAME --benchmark_repetitions=10 --benchmark_report_aggregates_only=true` |
| Capturing a baseline JSON for diffing | `--benchmark_format=json --benchmark_out=baseline.json --benchmark_repetitions=10` |
| Quick compile-link smoke test | `--benchmark_dry_run` |
| Investigating a high-CV benchmark | `taskset -c 0 ... --benchmark_repetitions=20 --benchmark_min_warmup_time=2.0` |
| Diffing two runs | `tools/compare.py benchmarks before.json after.json` |

