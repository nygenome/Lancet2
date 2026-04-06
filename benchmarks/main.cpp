#include "benchmark/benchmark.h"
#ifndef __APPLE__
#include "mimalloc-override.h"  // NOLINT(misc-include-cleaner)
#endif

BENCHMARK_MAIN();  // NOLINT
