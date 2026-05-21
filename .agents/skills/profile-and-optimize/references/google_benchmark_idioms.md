# google/benchmark idiom cheat sheet

This reference is loaded on demand by the `profile-and-optimize` skill when
writing or interpreting microbenchmarks. The SKILL itself stays slim;
this file holds the patterns that get misremembered or under-used.


Lancet2's microbenchmarks live in `benchmarks/` and use google/benchmark. The library docs are at [github.com/google/benchmark/tree/v1.9.5/docs](https://github.com/google/benchmark/tree/v1.9.5/docs); this section surfaces the idioms most often misremembered or under-used. Read it when you find yourself writing repetitive benchmark scaffolding that a built-in feature would replace, or when your benchmark numbers don't make sense.

### Anatomy of a benchmark

The basic shape is a function taking `benchmark::State&`, with the work in a range-for loop over the state. The library decides the iteration count automatically based on observed time per iteration:

```cpp
#include <benchmark/benchmark.h>

static void BM_HashKmer(benchmark::State& state) {
    auto kmer = MakeKmer(state.range(0));
    for (auto _ : state) {
        auto h = kmer.Hash();
        benchmark::DoNotOptimize(h);
    }
    state.SetBytesProcessed(int64_t(state.iterations()) * state.range(0));
}
BENCHMARK(BM_HashKmer)->Range(15, 31);
```

Two things to know about the loop. First, the work inside `for (auto _ : state)` is what gets timed. Setup before the loop is free; teardown after the loop is free. Second, the compiler will happily delete a loop whose result is unused — see "Preventing optimization" below.

### Passing arguments

A benchmark family is parameterized via the registration chain. Several patterns:

```cpp
BENCHMARK(BM_Foo)->Arg(8)->Arg(64)->Arg(512);              // explicit values
BENCHMARK(BM_Foo)->Range(8, 8<<10);                        // sparse range, /8 mult by default
BENCHMARK(BM_Foo)->RangeMultiplier(2)->Range(8, 8<<10);    // sparse range, /2 mult
BENCHMARK(BM_Foo)->DenseRange(0, 1024, 128);               // dense: 0,128,256,...,1024
BENCHMARK(BM_Foo)->Args({1024, 128})->Args({2048, 128});   // multi-dim explicit
BENCHMARK(BM_Foo)->Ranges({{1<<10, 8<<10}, {128, 512}});   // multi-dim sparse
BENCHMARK(BM_Foo)->ArgsProduct({{1024, 2048}, {20, 40}});  // cartesian product
BENCHMARK(BM_Foo)->Apply(CustomArguments);                 // programmatic
```

Inside the body, the values are accessed as `state.range(0)`, `state.range(1)`, etc. Names of generated benchmarks include the args (e.g., `BM_Foo/1024/128`), which makes them filterable.

For arbitrary non-integer arguments — strings, references to fixtures, anything that won't fit in `int64_t` — use `BENCHMARK_CAPTURE`:

```cpp
template <class... Args>
void BM_RunOn(benchmark::State& state, std::string_view fixture, Args&&... args) {
    auto data = LoadFixture(fixture);
    for (auto _ : state) {
        Process(data, args...);
    }
}
BENCHMARK_CAPTURE(BM_RunOn, na12878_chr1, "na12878_chr1.bam", 21);
BENCHMARK_CAPTURE(BM_RunOn, hcc1395_chr4, "hcc1395_chr4.bam", 21);
```

The `test_case_name` (`na12878_chr1`, `hcc1395_chr4` above) becomes part of the benchmark's display name.

### Templated benchmarks

Pass template parameters with the new C++17+ syntax:

```cpp
template <class Container>
void BM_Insert(benchmark::State& state) {
    Container c;
    for (auto _ : state) {
        c.insert(c.end(), state.range(0));
    }
}
BENCHMARK(BM_Insert<std::vector<int>>)->Arg(1000);
BENCHMARK(BM_Insert<std::list<int>>)->Arg(1000);
BENCHMARK(BM_Insert<std::deque<int>>)->Arg(1000);
```

The legacy `BENCHMARK_TEMPLATE(BM_Insert, std::vector<int>)` syntax also works; the inline form is preferred in new code.

### Fixtures

For setup/teardown that's too expensive to do per iteration but per benchmark is fine, derive from `benchmark::Fixture`:

```cpp
class GraphFixture : public benchmark::Fixture {
public:
    void SetUp(benchmark::State& state) override {
        // Built once per benchmark family member, before timing begins.
        mGraph = BuildTestGraph(state.range(0));
    }
    void TearDown(benchmark::State& state) override {}
protected:
    cbdg::Graph mGraph;
};

BENCHMARK_F(GraphFixture, BM_Traverse)(benchmark::State& state) {
    for (auto _ : state) {
        benchmark::DoNotOptimize(mGraph.Traverse());
    }
}
BENCHMARK_REGISTER_F(GraphFixture, BM_Traverse)->Arg(100)->Arg(1000)->Arg(10000);
```

Three macros for fixture benchmarks:
- `BENCHMARK_F(Class, Method)` — defines and registers in one shot.
- `BENCHMARK_DEFINE_F(Class, Method)` followed by `BENCHMARK_REGISTER_F(Class, Method)->...` — splits definition and registration so you can chain options on the registration.
- `BENCHMARK_TEMPLATE_F(Class, Method, T1, T2, ...)` — same for templated fixtures.

`SetUp` runs once before the timed loop, `TearDown` after. They are not part of the measurement.

### Pause/Resume timing inside the loop

Sometimes per-iteration setup is unavoidable — for example, generating a fresh random input each time. `state.PauseTiming()` and `state.ResumeTiming()` exclude code from the measurement:

```cpp
static void BM_SortRandom(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        auto data = GenerateRandomVec(state.range(0));
        state.ResumeTiming();
        std::sort(data.begin(), data.end());
    }
}
```

The library docs warn that `PauseTiming`/`ResumeTiming` has high overhead — it stops and starts the timer, which costs time that doesn't fit inside iterations. For benchmarks where the per-iteration setup is small, prefer `UseManualTime` (below) or batched setup before the loop.

### Manual timing

For work where neither CPU time nor wall clock is the right measurement, take the timing yourself and report it:

```cpp
static void BM_HardwareSpecific(benchmark::State& state) {
    for (auto _ : state) {
        auto t0 = std::chrono::high_resolution_clock::now();
        DoTheThing();
        auto t1 = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration<double>(t1 - t0).count();
        state.SetIterationTime(elapsed);
    }
}
BENCHMARK(BM_HardwareSpecific)->UseManualTime();
```

`UseManualTime()` tells the library to ignore its own timer and use what `SetIterationTime` reports. This is the right tool when measuring something that happens partly outside CPU time (GPU calls, I/O on a different device, sleeps).

### Preventing optimization

The single most common benchmark bug: the optimizer deletes the benchmarked code because nothing depends on its result. Two helpers:

- **`benchmark::DoNotOptimize(expr)`** — forces the *result* of `expr` to be stored in a register or memory, so the compiler can't dead-code-eliminate it. This does NOT prevent constant folding of the expression; if `expr` is computable at compile time, the compiler may still pre-compute it. Pass an actual variable or a non-trivial expression.

- **`benchmark::ClobberMemory()`** — forces the compiler to flush all pending writes to memory. Use after a write that the compiler might otherwise consider unnecessary.

The canonical pattern for benchmarking a write into a container:

```cpp
static void BM_VectorPushBack(benchmark::State& state) {
    for (auto _ : state) {
        std::vector<int> v;
        v.reserve(1);
        auto* data = v.data();
        benchmark::DoNotOptimize(data);  // escape data so its writes are observable
        v.push_back(42);
        benchmark::ClobberMemory();      // force the write to actually happen
    }
}
```

If your benchmark numbers look suspiciously fast, suspect missing `DoNotOptimize`/`ClobberMemory` first — well before "the algorithm is just that fast."

### Custom counters

Beyond the default time/iteration columns, the library exposes a per-benchmark counter map:

```cpp
static void BM_ProcessReads(benchmark::State& state) {
    int64_t reads = 0;
    for (auto _ : state) {
        reads += ProcessBatch(state.range(0));
    }
    state.counters["reads/sec"] =
        benchmark::Counter(reads, benchmark::Counter::kIsRate);
    state.counters["reads"] = reads;  // total, no flags
}
```

Counter flags (combinable with `|`):
- `kIsRate` — divide by total time, present as throughput.
- `kAvgThreads` — divide by thread count.
- `kAvgThreadsRate` — both above.
- `kIsIterationInvariantRate` — multiply by iteration count, then divide by time.
- `kInvert` — invert (e.g., turn rate into seconds-per-thing).

For Lancet2, `bytes_per_second` (via `state.SetBytesProcessed`) and `items_per_second` (via `state.SetItemsProcessed`) are the conventional throughput metrics; reach for custom counters when you need something else.

### Setup/Teardown callbacks

Distinct from fixture `SetUp`/`TearDown` — these are global-like callbacks invoked once per registered family member:

```cpp
static void DoSetup(const benchmark::State& state) { /* ... */ }
static void DoTeardown(const benchmark::State& state) { /* ... */ }

BENCHMARK(BM_func)
    ->Arg(1)->Arg(3)
    ->Threads(16)->Threads(32)
    ->Setup(DoSetup)->Teardown(DoTeardown);
```

The callbacks run once per (arg, thread-count) combination — in the example above, four invocations.

### Disabling a benchmark

Prefix the function name with `DISABLED_`:

```cpp
static void DISABLED_BM_Slow(benchmark::State& state) { /* ... */ }
BENCHMARK(DISABLED_BM_Slow);
```

Disabled benchmarks register but skip at runtime. The display name keeps the prefix, so it's visually obvious which ones are off.

