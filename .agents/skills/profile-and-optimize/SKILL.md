---
name: profile-and-optimize
description: Use when investigating Lancet2 performance — slow runs, hot functions, regressions, or candidate optimizations — OR when figuring out how to write, run, or interpret a microbenchmark or CPU profile. Trigger on "profile", "optimize", "slow", "hotspot", "performance", "wall-clock", "benchmark", "DoNotOptimize", "ClobberMemory", "BENCHMARK_F", "benchmark filter", "benchmark variance", "interpret this profile", "what does pprof show", "pprof focus", "analyze_profile", "compare profiles". Builds the profiling tree, runs against real data, analyzes via the project's `analyze_profile.py` (a pprof wrapper), optionally delegates to the perf-analyst subagent, validates each change with both pipeline wall-clock and the google/benchmark microbenchmark suite under `benchmarks/`. The procedural detail — google/benchmark idioms, bench-binary invocation, gperftools internals, pprof and the analyze_profile.py wrapper — lives in `references/` files loaded on demand. The SKILL is the trigger and the workflow; the references hold the depth.
allowed-tools: Read, Glob, Grep, Edit, Write, Bash
---

# Profile-and-optimize on Lancet2

This skill formalizes the full performance-improvement workflow. It composes with the `perf-analyst` subagent (for analysis when context isolation matters) and includes the canonical benchmark procedure for unit-level measurement. The discipline is to never accept "this should be faster" without a benchmark number; performance intuition is wrong often enough that unverified optimizations produce more regressions than gains.

The procedural detail — google/benchmark idioms, bench-binary command-line surface, gperftools internals, pprof and the project's `analyze_profile.py` wrapper — lives in four reference files under `references/`. Each step below names the reference to consult when you need depth. The SKILL stays slim so it's cheap to load on every performance-tagged trigger; the references load only when the work demands them.

## Step 1 — Build the profiling tree

The project ships a complete profiling pipeline. Configure and build the profiling tree, which uses `RelWithDebInfo` for source-line annotation and links gperftools statically:

```bash
pixi run configure-profile
pixi run build-profile
```

The resulting binary is `cmake-build-relwithdebinfo/Lancet2`. The `LANCET_PROFILE_MODE` build option enables CPU profiling; `platform_checks.cmake` enforces the build type for you.

## Step 2 — Capture a profile against real data

Run the binary on real data, against a region representative of your typical workload. **Lancet2 does not honor the `CPUPROFILE` environment variable** — the profile filename and the sampling settings (frequency 250 Hz, per-thread timers) are baked into the source at `src/lancet/cli/pipeline_runner.cpp`. The `setenv` calls there fire before `ProfilerStart` and are the only way these settings take effect, because gperftools reads its env vars during static initialization. To change the profile filename or the sampling frequency, edit the source and rebuild — see `references/gperftools_profiler.md` for the full picture.

Practical consequence: pick the working directory carefully before invoking the binary, because that's where `Lancet.cpu_profile.<TIMESTAMP>.bin` lands.

```bash
mkdir -p /tmp/lancet2-profile && cd /tmp/lancet2-profile
"$OLDPWD/cmake-build-relwithdebinfo/Lancet2" pipeline \
  --tumor $LANCET_TEST_SOMATIC_TUMOR \
  --normal $LANCET_TEST_SOMATIC_NORMAL \
  --reference $LANCET_TEST_SOMATIC_REFERENCE \
  --region $LANCET_TEST_SOMATIC_REGION \
  --num-threads $(nproc) \
  --out-vcfgz /tmp/lancet2-profile/profile_run.vcf.gz
cd "$OLDPWD"
LANCET_PROFILE_BIN=$(ls -t /tmp/lancet2-profile/Lancet.cpu_profile.*.bin | head -1)
echo "Profile: $LANCET_PROFILE_BIN"
```

Use `$LANCET_TEST_SOMATIC_REGION` (the full chr4 sweep) rather than `$LANCET_TEST_SOMATIC_REGION_SMALL`; profiling is noisy on tiny inputs, and a region with at least a few hundred variants gives a stable profile.

If the run is too short to produce a useful profile, increase the region size or chain multiple regions. If it is too long, reduce the thread count to one for a deterministic single-threaded profile, then capture a multi-threaded one separately.

## Step 3 — Analyze with the project's pprof wrapper

The project has a sophisticated profile-analysis script. Run it in the profiling environment:

```bash
pixi run -e profiling analyze-profile "$LANCET_PROFILE_BIN" -- --view overview
pixi run -e profiling analyze-profile "$LANCET_PROFILE_BIN" -- --view top --top 50
pixi run -e profiling analyze-profile "$LANCET_PROFILE_BIN" -- --view tree --focus "lancet::caller"
pixi run -e profiling analyze-profile "$LANCET_PROFILE_BIN" -- --view hotpaths
```

The `overview` view tells you total samples and top namespaces. The `top` view identifies individual functions by self time. The `tree` view shows call-graph structure under a chosen subtree. The `hotpaths` view ranks the hottest call chains end-to-end. For visual exploration, generate an HTML report:

```bash
pixi run -e profiling analyze-profile "$LANCET_PROFILE_BIN" -- --html /tmp/profile.html
```

## Step 4 — Hand the analysis to perf-analyst

Capture the relevant excerpts from the views above and invoke the subagent:

```
Use the perf-analyst subagent to interpret this profile and propose optimizations:

[paste --view overview output]

[paste --view top --top 30 output]

[paste --view tree --focus output for the highest-cost namespace]
```

The subagent will return three to five candidates with explicit trade-off analysis, the smallest test that would validate each, and a recommended order based on cost-to-implement and confidence in the gain.

## Step 5 — Save the baseline summary

Before implementing any optimization, save a summary of the baseline so you can diff against it later:

```bash
pixi run -e profiling analyze-profile "$LANCET_PROFILE_BIN" -- \
  --binary cmake-build-relwithdebinfo/Lancet2 \
  --save-summary baseline
```

This stores the summary keyed as "baseline" for later cross-binary diff. The `--save-summary` flow is the canonical optimization workflow because it survives binary rebuilds (a regular `--diff-base` comparison fails when the binary changes).

## Step 6 — Implement ONE candidate

Pick the highest-confidence, lowest-cost candidate from the perf-analyst recommendations. Use plan mode in the calling session, with Opus, and work in a worktree if the change is non-trivial:

```bash
git worktree add ../lancet2-perf-candidate1 perf-candidate1
cd ../lancet2-perf-candidate1
```

Implement the change in `src/lancet/<layer>/`, respecting the layer-direction rule. Run `pixi run test` to confirm correctness; correctness must come before performance. Run `pixi run lint-check` to confirm convention adherence.

## Step 7 — Re-profile and diff

Rebuild the profiling tree and capture a new profile. As in step 2, the filename and sampling settings are baked into the binary; the only knob the user has at runtime is the working directory:

```bash
pixi run build-profile
mkdir -p /tmp/lancet2-after && cd /tmp/lancet2-after
"$OLDPWD/cmake-build-relwithdebinfo/Lancet2" pipeline \
  --tumor $LANCET_TEST_SOMATIC_TUMOR \
  --normal $LANCET_TEST_SOMATIC_NORMAL \
  --reference $LANCET_TEST_SOMATIC_REFERENCE \
  --region $LANCET_TEST_SOMATIC_REGION \
  --num-threads $(nproc) \
  --out-vcfgz /tmp/lancet2-after/after_run.vcf.gz
cd "$OLDPWD"
LANCET_AFTER_BIN=$(ls -t /tmp/lancet2-after/Lancet.cpu_profile.*.bin | head -1)
echo "After profile: $LANCET_AFTER_BIN"
```

Save the new summary and compare against baseline:

```bash
pixi run -e profiling analyze-profile "$LANCET_AFTER_BIN" -- \
  --binary cmake-build-relwithdebinfo/Lancet2 \
  --save-summary candidate1
pixi run -e profiling analyze-profile -- --diff-tag baseline candidate1
```

The diff output identifies which functions got faster, which got slower (always check for regressions), and the net wall-clock impact.

## Step 8 — Validate VCF output identity

A performance change must not change the variant calls. Compare the two VCFs:

```bash
pixi run -e hts-tools bcftools view /tmp/lancet2-profile/profile_run.vcf.gz | grep -v '^##' > /tmp/before.vcf
pixi run -e hts-tools bcftools view /tmp/lancet2-after/after_run.vcf.gz   | grep -v '^##' > /tmp/after.vcf
diff /tmp/before.vcf /tmp/after.vcf
```

For correctness-critical changes (anything touching scoring math, graph traversal, or genotype likelihoods), the diff must be empty. For changes that only affect throughput (memory layout, allocation patterns, threading), the diff is also expected to be empty for the same input. If the diff is not empty, the change has changed semantics and must be reverted or justified explicitly.

## Step 9 — Cross-validate with the benchmark suite

Beyond the pipeline-level wall clock, the project has a google/benchmark suite under `benchmarks/` that gives microsecond-resolution measurement with explicit statistical context. The pipeline-level wall clock is the right metric for end-user impact, but it is too coarse and noisy to validate small targeted changes; the benchmark suite is the right metric for unit-level confidence.

**When to benchmark vs. when the pipeline measurement is enough.** Use the benchmark suite when the change targets a specific function or class (e.g., a graph-builder hotspot or a scoring function) and you want microsecond-resolution before/after numbers. Stick with the pipeline measurement (the `--diff-tag` step above) when the change is broad enough that wall clock is the right metric, or when the function does not have a benchmark and writing one is more work than it is worth.

**Locating or writing the benchmark.** Benchmarks live in `benchmarks/` and mirror the structure of `src/lancet/`. Glob the directory first; extending an existing benchmark with a new `BENCHMARK` block is preferable to creating a new file. If no benchmark exists, copy the closest existing one and adapt it — read at least two existing benchmarks before writing yours so you pick up project conventions on input fixtures, randomization to avoid CPU-cache hot-spotting, and benchmark labeling.

**Building and running.** The benchmark executable is built when `LANCET_BENCHMARKS=ON` is passed to CMake. The pixi `bench` task handles this. RelWithDebInfo with `-DLANCET_BENCHMARKS=ON` is the right configuration for measurement work — optimizer-on performance with enough debug information for `perf` correlation.

```bash
pixi run --quiet bench
```

**Baseline and after.** Capture the baseline on `main` (or the change's branch base) before applying the change, then capture the after on the change branch:

```bash
git switch main
pixi run --quiet bench > /tmp/bench_baseline.txt 2>&1
git switch your-change-branch
pixi run --quiet bench > /tmp/bench_after.txt 2>&1
```

Run each side at least twice and confirm consistency. Single-run comparisons are unreliable because of frequency scaling, thermal throttling, and noisy neighbors.

**Interpreting the delta.** google/benchmark's report includes per-benchmark mean, median, standard deviation, and coefficient of variation when run with `--benchmark_repetitions=N`. A meaningful improvement is one where the after-mean is more than two standard deviations below the before-mean and the relative change exceeds the noise floor. If the improvement is below five percent, treat it as noise unless many runs are averaged. If it is above twenty percent, double-check that the change actually applied (a common error is benchmarking the wrong build directory). If a different benchmark unexpectedly regressed, investigate that before celebrating the targeted improvement. If the coefficient of variation on any benchmark exceeds about ten percent, the measurement environment is too noisy — pin the CPU governor, isolate cores, or increase iterations. The library's `tools/compare.py` can run a Mann-Whitney U-test on two JSON outputs and report deltas with significance; see `references/benchmark_invocation.md` (the "Comparing runs with `compare.py`" subsection) for the invocation and the recipe table.

**When the benchmark suite does not cover what you need.** If the function being changed has no benchmark and no existing benchmark is close enough to adapt, write a new benchmark file before making the performance change. The benchmark commit can land separately on `main` ahead of the optimization, which gives you a stable baseline and benefits future optimizers in the same area.

## Step 10 — Commit with measurement evidence

Use the conventional-commit prefix `perf:` and include the measured improvement in the message:

```
perf: hoist sample mask to outer loop in graph builder

Reduces inner-loop allocations from O(n*k) to O(n).
Wall clock on `<region>` (10x normal/tumor): 142s → 128s (-9.9%).
Benchmark "graph_build_complex": 1.84ms/iter → 1.62ms/iter (-12%).
VCF output identical (bcftools diff empty).
```

The CHANGELOG entry generated by git-chglog will surface this for future readers, who can decide whether the change is worth keeping if a future regression points back to it.

## Step 11 — Move to the next candidate or stop

If there is an obvious second candidate from the perf-analyst output, repeat from step 6 with the next worktree. If not, stop. Diminishing returns set in quickly; three optimizations in a row that each save 10 percent are far better than ten that each save 1 percent.

## When NOT to use this skill

Do not use this skill for micro-optimization (function inlining, loop unrolling, individual instruction-level changes); the compiler is better at these and the noise floor of measurement makes the gain unverifiable. Do not use it for changes whose purpose is not performance (refactoring, naming, modernization); those go through the standard add-cpp-test workflow. Do not use it without representative input data; profiling on a tiny region produces conclusions that do not survive contact with real workloads.

## References

These are loaded on demand when a step in the workflow references them.

- `references/google_benchmark_idioms.md` — `BENCHMARK_F`, `DoNotOptimize`, `ClobberMemory`, `Range`/`RangeMultiplier`/`DenseRange`, fixtures, pause/resume vs manual timing, custom counters, setup/teardown callbacks, disabling. Read this when writing a new benchmark or extending an existing one (Step 9).
- `references/benchmark_invocation.md` — every `--benchmark_*` flag and its env-var equivalent, subset selection, repetitions and statistics, output formats and JSON capture for `compare.py`, dry-run and list mode, variance-reduction patterns (governor pinning, `taskset`, warmup), and a recipe table for common scenarios. Read this when running the benchmark binary directly (Step 9).
- `references/gperftools_profiler.md` — how Lancet2 starts/stops profiling (settings baked in `pipeline_runner.cpp`), what gets sampled, wall-time vs CPU-time, upstream caveats, and Lancet2-specific pitfalls (forgetting the profiling build, looking in the wrong directory, trying to set `CPUPROFILE` from the shell). Read this when modifying profiling settings or interpreting profile capture (Steps 2 and 7).
- `references/pprof_and_analyze.md` — the wrapper's view system (`overview`, `top`, `tree`, `hotpaths`, `compare`, `--save-summary`, `--diff-tag`, HTML rendering), raw pprof reference, interactive shell commands, and the rule for when to use the wrapper vs raw pprof. Read this when analyzing profiles (Steps 3 and 7).
