// ============================================================================
// Repeat Benchmark — Hamming distance and repeat detection performance
//
// Part 1: Hamming distance implementation comparison
//   Benchmarks six approaches to byte-level Hamming distance:
//     1. HammingDistIntrinsics (production) — explicit SIMD from repeat.h
//     2. HammingDistAutoVec   — u8 batch accumulation (compiler auto-vectorized)
//     3. HammingDistUsizeAccum — usize accumulator (widening bottleneck)
//     4. HammingDistU32Accum  — uint32_t accumulator (partial fix)
//     5. HammingDistSWAR      — manual SWAR with popcount
//
// Part 2: HasRepeat with varying entropy and mismatch thresholds
//   Tests two implementations:
//     1. HasRepeatIntrinsics (production) — IsWithinHammingDist + early exit
//     2. HasRepeatAutoVec                — HammingDistAutoVec + full distance
//   Across three sequence entropy levels:
//     Low    — homopolymer-rich (many repeats, fast short-circuit)
//     Medium — random DNA (realistic k-mer distribution)
//     High   — all 4 bases uniformly shuffled (few repeats, worst case)
// ============================================================================

#include "lancet/base/repeat.h"

#include "lancet/base/assert.h"
#include "lancet/base/sliding.h"
#include "lancet/base/types.h"

#include "absl/container/flat_hash_set.h"
#include "absl/random/distributions.h"
#include "absl/types/span.h"
#include "benchmark/benchmark.h"

#include <algorithm>
#include <array>
#include <bit>
#include <functional>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

// ── Weighted index pick over a cumulative-threshold table ──────────────────
// Stand-in for `std::discrete_distribution` (which has no 1:1 absl
// replacement). Encapsulates the half-open `absl::Uniform<u32>(gen, 0,
// cumulative_weights.back())` invariant so callers cannot accidentally
// break loop termination by switching to `absl::IntervalClosed`. Picks
// an index in `[0, N)` proportional to the per-bucket weights implied
// by the cumulative-threshold table; the trailing `return N - 1`
// handles the (statistically-impossible-with-half-open) tie at the
// upper bound, keeping the helper safe even if the draw form is later
// changed by a maintainer.
template <usize n>
[[nodiscard]] inline auto WeightedPick(std::array<u32, n> const& cumulative_weights,
                                       std::mt19937_64& gen) -> usize {
  auto const draw = absl::Uniform<u32>(gen, 0, cumulative_weights.back());
  for (usize idx = 0; idx < n - 1; ++idx) {
    if (draw < cumulative_weights[idx]) return idx;
  }
  return n - 1;
}

// ── Sequence generators with controllable entropy ──────────────────────────

// Low entropy: heavily biased toward A (75% A, ~8% each C/G/T).
// Produces homopolymer runs and frequent exact/approximate repeats.
[[nodiscard]] inline auto GenerateLowEntropySequence(usize const seq_len) -> std::string {
  static constexpr std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};

  // Fixed seed for reproducible benchmarks across runs
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937_64 generator(42);

  // Heavily biased toward A — 75% A, ~8% each for C/G/T. Cumulative
  // boundaries derived from the {75, 8, 8, 9} weight vector.
  static constexpr std::array<u32, 4> CUM_WEIGHTS = {75, 83, 91, 100};
  std::string result(seq_len, 'N');

  for (usize idx = 0; idx < seq_len; ++idx) {
    result[idx] = BASES.at(WeightedPick(CUM_WEIGHTS, generator));
  }

  return result;
}

// Medium entropy: uniform random DNA — the baseline for typical genomic sequence.
[[nodiscard]] inline auto GenerateMediumEntropySequence(usize const seq_len) -> std::string {
  static constexpr std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};

  // Fixed seed for reproducible benchmarks across runs
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937_64 generator(137);

  std::string result(seq_len, 'N');

  for (usize idx = 0; idx < seq_len; ++idx) {
    result[idx] = BASES.at(absl::Uniform<usize>(absl::IntervalClosed, generator, 0, 3));
  }

  return result;
}

// High entropy: maximally diverse — all 16 dinucleotides appear with equal
// probability by cycling a shuffled alphabet.  Minimizes k-mer collisions.
[[nodiscard]] inline auto GenerateHighEntropySequence(usize const seq_len) -> std::string {
  static constexpr std::array<char, 16> CYCLE = {'A', 'C', 'G', 'T', 'C', 'A', 'T', 'G',
                                                 'G', 'T', 'A', 'C', 'T', 'G', 'C', 'A'};

  // Fixed seed for reproducible benchmarks across runs
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937_64 generator(271);

  std::string result(seq_len, 'N');

  for (usize idx = 0; idx < seq_len; ++idx) {
    result[idx] = CYCLE.at(absl::Uniform<usize>(absl::IntervalClosed, generator, 0, 15));
  }

  return result;
}

// ════════════════════════════════════════════════════════════════════════════
// Baseline implementations for benchmarking (not used in production)
// ════════════════════════════════════════════════════════════════════════════

// ── Auto-vectorized HammingDist (former production) ────────────────────────
// u8 batch accumulation proves to the compiler's value-range analysis that no
// byte lane overflows, allowing the auto-vectorizer to emit vpcmpeqb + vpsadbw
// on AVX2 without widening.  4.7× fewer SIMD instructions than usize accumulator.
// Replaced in production by explicit SIMD intrinsics with overlapping-tail loads.
auto HammingDistAutoVec(std::string_view first, std::string_view second) -> usize {
  LANCET_ASSERT(first.length() == second.length())

  usize result = 0;
  auto const length = first.length();
  usize idx = 0;

  // 255 is the largest count a u8 can hold without overflow (each iteration
  // adds 0 or 1).  One scalar widen per 255 bytes is negligible overhead.
  while (idx < length) {
    u8 batch_sum = 0;
    auto const batch_end = std::min(idx + 255, length);
    for (; idx < batch_end; ++idx) {
      batch_sum += static_cast<u8>(first[idx] != second[idx]);
    }
    result += batch_sum;
  }

  return result;
}

// ── HasRepeat using auto-vectorized HammingDist (former production) ────────
// Full-distance comparison: always computes the entire Hamming distance before
// comparing against the threshold.  No early exit within the distance function.
// Replaced in production by IsWithinHammingDist which early-exits per SIMD chunk.
auto HasRepeatAutoVec(absl::Span<std::string_view const> kmers, usize const max_mismatches)
    -> bool {
  if (max_mismatches == 0) {
    absl::flat_hash_set<std::string_view> seen;
    seen.reserve(kmers.size());
    for (auto const kmer : kmers) {
      if (auto [iter, inserted] = seen.insert(kmer); !inserted) return true;
    }
    return false;
  }

  auto const num_kmers = kmers.size();
  for (usize i = 0; i < num_kmers; ++i) {
    for (usize j = i + 1; j < num_kmers; ++j) {
      if (HammingDistAutoVec(kmers[i], kmers[j]) <= max_mismatches) return true;
    }
  }
  return false;
}

// ── Rejected approach 1: usize accumulator ─────────────────────────────────
// The widening bottleneck: each 1-byte comparison result is zero-extended to
// 8 bytes before accumulation.  GCC 15.2 emits 46 vpmovzx instructions and
// 168 total SIMD instructions per vectorized iteration on AVX2.
auto HammingDistUsizeAccum(std::string_view first, std::string_view second) -> usize {
  usize result = 0;
  auto const length = first.length();
  for (usize idx = 0; idx < length; ++idx) {
    result += static_cast<usize>(first[idx] != second[idx]);
  }
  return result;
}

// ── Rejected approach 2: uint32_t accumulator ──────────────────────────────
// Halves the widening cost (vpmovzxbd packs 8 elements vs vpmovzxbq's 4),
// but still leaves 18 vpmovzx instructions on GCC.  No effect on Clang.
auto HammingDistU32Accum(std::string_view first, std::string_view second) -> usize {
  u32 local_res = 0;
  auto const length = first.length();
  for (usize idx = 0; idx < length; ++idx) {
    local_res += static_cast<u32>(first[idx] != second[idx]);
  }
  return static_cast<usize>(local_res);
}

// ── Rejected approach 3: manual SWAR with popcount ─────────────────────────
// Processes 8 bytes per iteration using shift+OR+mask to fold each byte lane
// into a single bit, then popcount.  Defeats the auto-vectorizer entirely:
// the serial dependency chain and scalar popcount prevent SIMD widening.
// Also contains three instances of undefined behavior (buffer over-read,
// alignment violation, strict aliasing violation via reinterpret_cast).
// Kept here purely as a benchmark baseline — do NOT use in production.
auto HammingDistSWAR(std::string_view first, std::string_view second) -> usize {
  usize result = 0;

  auto const num_words = (first.length() >> 3);
  auto const rem_words = static_cast<unsigned long long>(first.length() & 7);

  // SWAR (SIMD Within A Register) reinterprets the byte buffer as a packed 64-bit word stream;
  // the cast is the conventional idiom for SWAR primitives.
  // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* aptr = reinterpret_cast<unsigned long long const*>(first.data());
  auto const* bptr = reinterpret_cast<unsigned long long const*>(second.data());
  // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

  for (usize idx = 0; idx < num_words; ++idx) {
    auto val = (aptr[idx] ^ bptr[idx]);
    val |= val >> 4;
    val &= 0x0f'0f'0f'0f'0f'0f'0f'0fULL;
    val |= val >> 2;
    val &= 0x33'33'33'33'33'33'33'33ULL;
    val |= val >> 1;
    val &= 0x55'55'55'55'55'55'55'55ULL;
    result += std::popcount(val);
  }

  if (rem_words > 0) {
    auto val = (aptr[num_words] ^ bptr[num_words]);
    val |= val >> 4;
    val &= 0x0f'0f'0f'0f'0f'0f'0f'0fULL;
    val |= val >> 2;
    val &= 0x33'33'33'33'33'33'33'33ULL;
    val |= val >> 1;
    val &= 0x55'55'55'55'55'55'55'55ULL;
    result += std::popcount((val & ((1ULL << (rem_words << 3ULL)) - 1ULL)));
  }

  return result;
}

// ── Benchmark infrastructure ───────────────────────────────────────────────

using GeneratorFn = std::string (*)(usize);

// ════════════════════════════════════════════════════════════════════════════
// Part 1: Hamming distance — implementation comparison
// ════════════════════════════════════════════════════════════════════════════

// ── Production: explicit SIMD intrinsics (AVX2/NEON) ───────────────────────
template <GeneratorFn generator>
void BenchHammingDistIntrinsics(benchmark::State& state) {
  auto const seq_len = static_cast<usize>(state.range(0));
  std::string const first = generator(seq_len);
  std::string const second = generator(seq_len);

  // google/benchmark idiom: `_` is the conventional name for the unused iteration variable.
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    auto result = lancet::base::HammingDist(first, second);
    benchmark::DoNotOptimize(result);
  }
}

// ── Former production: auto-vectorized u8 batch accumulation ───────────────
template <GeneratorFn generator>
void BenchHammingDistAutoVec(benchmark::State& state) {
  auto const seq_len = static_cast<usize>(state.range(0));
  std::string const first = generator(seq_len);
  std::string const second = generator(seq_len);

  // google/benchmark idiom: `_` is the conventional name for the unused iteration variable.
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    auto result = HammingDistAutoVec(first, second);
    benchmark::DoNotOptimize(result);
  }
}

// ── Rejected: usize accumulator ────────────────────────────────────────────
template <GeneratorFn generator>
void BenchHammingDistUsizeAccum(benchmark::State& state) {
  auto const seq_len = static_cast<usize>(state.range(0));
  std::string const first = generator(seq_len);
  std::string const second = generator(seq_len);

  // google/benchmark idiom: `_` is the conventional name for the unused iteration variable.
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    auto result = HammingDistUsizeAccum(first, second);
    benchmark::DoNotOptimize(result);
  }
}

// ── Rejected: u32 accumulator ──────────────────────────────────────────────
template <GeneratorFn generator>
void BenchHammingDistU32Accum(benchmark::State& state) {
  auto const seq_len = static_cast<usize>(state.range(0));
  std::string const first = generator(seq_len);
  std::string const second = generator(seq_len);

  // google/benchmark idiom: `_` is the conventional name for the unused iteration variable.
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    auto result = HammingDistU32Accum(first, second);
    benchmark::DoNotOptimize(result);
  }
}

// ── Rejected: SWAR + popcount ──────────────────────────────────────────────
template <GeneratorFn generator>
void BenchHammingDistSWAR(benchmark::State& state) {
  auto const seq_len = static_cast<usize>(state.range(0));
  std::string const first = generator(seq_len);
  std::string const second = generator(seq_len);

  // google/benchmark idiom: `_` is the conventional name for the unused iteration variable.
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    auto result = HammingDistSWAR(first, second);
    benchmark::DoNotOptimize(result);
  }
}

// ════════════════════════════════════════════════════════════════════════════
// Part 2: HasRepeat — entropy × mismatch threshold matrix
//
// state.range(0) = sequence length (bp)
// state.range(1) = k-mer window size
// state.range(2) = max_mismatches (0 = exact repeat, 1–8 = approximate)
// ════════════════════════════════════════════════════════════════════════════

// ── Production: SIMD intrinsics with early-exit ────────────────────────────
template <GeneratorFn generator>
void BenchHasRepeatIntrinsics(benchmark::State& state) {
  auto const seq_len = static_cast<usize>(state.range(0));
  auto const window = static_cast<usize>(state.range(1));
  auto const max_mismatches = static_cast<usize>(state.range(2));
  std::string const seq = generator(seq_len);
  auto const kmers = lancet::base::SlidingView(seq, window);

  // google/benchmark idiom: `_` is the conventional name for the unused iteration variable.
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    auto result = lancet::base::HasRepeat(absl::MakeConstSpan(kmers), max_mismatches);
    benchmark::DoNotOptimize(result);
  }

  state.SetLabel(std::to_string(max_mismatches) + " mismatches");
}

// ── Former production: auto-vectorized full-distance HasRepeat ──────────────
template <GeneratorFn generator>
void BenchHasRepeatAutoVec(benchmark::State& state) {
  auto const seq_len = static_cast<usize>(state.range(0));
  auto const window = static_cast<usize>(state.range(1));
  auto const max_mismatches = static_cast<usize>(state.range(2));
  std::string const seq = generator(seq_len);
  auto const kmers = lancet::base::SlidingView(seq, window);

  // google/benchmark idiom: `_` is the conventional name for the unused iteration variable.
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    auto result = HasRepeatAutoVec(absl::MakeConstSpan(kmers), max_mismatches);
    benchmark::DoNotOptimize(result);
  }

  state.SetLabel(std::to_string(max_mismatches) + " mismatches");
}

// Parameterize: (seq_length, kmer_window, max_mismatches)
// seq_length=600 approximates a typical Lancet2 genomic window.
// kmer_window values from graph_params.h: min=13, max=127, step=6.
void ApplyHasRepeatArgs(benchmark::Benchmark* bench) {
  static constexpr std::array<i64, 4> KMER_SIZES = {13, 55, 97, 127};
  static constexpr i64 SEQ_LENGTH = 600;

  for (auto const kmer : KMER_SIZES) {
    // exact repeat (max_mismatches = 0)
    bench->Args({SEQ_LENGTH, kmer, 0});
    // approximate repeats: 1 through 8 mismatches
    for (i64 mm = 1; mm <= 8; ++mm) {
      bench->Args({SEQ_LENGTH, kmer, mm});
    }
  }
}

}  // namespace

// google/benchmark BENCHMARK() macros instantiate static registrar objects at namespace scope
// (cert-err58-cpp), allocate fluent-API state via raw new (owning-memory), are required to live
// at namespace scope outside an anonymous namespace so the macro emits external linkage symbols
// (use-anonymous-namespace), and use the library-defined short macro name (identifier-length).
// NOLINTBEGIN(cert-err58-cpp, cppcoreguidelines-owning-memory, readability-identifier-length, misc-use-anonymous-namespace)

// ── Hamming distance: k-mer sized inputs (11–121 bp), all entropy levels ──

// Production: explicit SIMD intrinsics
BENCHMARK(BenchHammingDistIntrinsics<GenerateLowEntropySequence>)->DenseRange(11, 121, 10);
BENCHMARK(BenchHammingDistIntrinsics<GenerateMediumEntropySequence>)->DenseRange(11, 121, 10);
BENCHMARK(BenchHammingDistIntrinsics<GenerateHighEntropySequence>)->DenseRange(11, 121, 10);

// Former production: auto-vectorized u8 batch
BENCHMARK(BenchHammingDistAutoVec<GenerateLowEntropySequence>)->DenseRange(11, 121, 10);
BENCHMARK(BenchHammingDistAutoVec<GenerateMediumEntropySequence>)->DenseRange(11, 121, 10);
BENCHMARK(BenchHammingDistAutoVec<GenerateHighEntropySequence>)->DenseRange(11, 121, 10);

// Rejected baselines
BENCHMARK(BenchHammingDistUsizeAccum<GenerateLowEntropySequence>)->DenseRange(11, 121, 10);
BENCHMARK(BenchHammingDistUsizeAccum<GenerateMediumEntropySequence>)->DenseRange(11, 121, 10);
BENCHMARK(BenchHammingDistUsizeAccum<GenerateHighEntropySequence>)->DenseRange(11, 121, 10);

BENCHMARK(BenchHammingDistU32Accum<GenerateLowEntropySequence>)->DenseRange(11, 121, 10);
BENCHMARK(BenchHammingDistU32Accum<GenerateMediumEntropySequence>)->DenseRange(11, 121, 10);
BENCHMARK(BenchHammingDistU32Accum<GenerateHighEntropySequence>)->DenseRange(11, 121, 10);

BENCHMARK(BenchHammingDistSWAR<GenerateLowEntropySequence>)->DenseRange(11, 121, 10);
BENCHMARK(BenchHammingDistSWAR<GenerateMediumEntropySequence>)->DenseRange(11, 121, 10);
BENCHMARK(BenchHammingDistSWAR<GenerateHighEntropySequence>)->DenseRange(11, 121, 10);

// ── Hamming distance: power-of-two stress test (8–2048 bytes) ─────────────

// Production: explicit SIMD intrinsics
BENCHMARK(BenchHammingDistIntrinsics<GenerateLowEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);
BENCHMARK(BenchHammingDistIntrinsics<GenerateMediumEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);
BENCHMARK(BenchHammingDistIntrinsics<GenerateHighEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);

// Former production: auto-vectorized u8 batch
BENCHMARK(BenchHammingDistAutoVec<GenerateLowEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);
BENCHMARK(BenchHammingDistAutoVec<GenerateMediumEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);
BENCHMARK(BenchHammingDistAutoVec<GenerateHighEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);

// Rejected baselines
BENCHMARK(BenchHammingDistUsizeAccum<GenerateLowEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);
BENCHMARK(BenchHammingDistUsizeAccum<GenerateMediumEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);
BENCHMARK(BenchHammingDistUsizeAccum<GenerateHighEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);

BENCHMARK(BenchHammingDistU32Accum<GenerateLowEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);
BENCHMARK(BenchHammingDistU32Accum<GenerateMediumEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);
BENCHMARK(BenchHammingDistU32Accum<GenerateHighEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);

BENCHMARK(BenchHammingDistSWAR<GenerateLowEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);
BENCHMARK(BenchHammingDistSWAR<GenerateMediumEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);
BENCHMARK(BenchHammingDistSWAR<GenerateHighEntropySequence>)
    ->RangeMultiplier(2)
    ->Range(2 << 2, 2 << 10);

// ── HasRepeat: production SIMD intrinsics with early exit ─────────────────
BENCHMARK(BenchHasRepeatIntrinsics<GenerateLowEntropySequence>)->Apply(ApplyHasRepeatArgs);
BENCHMARK(BenchHasRepeatIntrinsics<GenerateMediumEntropySequence>)->Apply(ApplyHasRepeatArgs);
BENCHMARK(BenchHasRepeatIntrinsics<GenerateHighEntropySequence>)->Apply(ApplyHasRepeatArgs);

// ── HasRepeat: former auto-vectorized full-distance comparison ────────────
BENCHMARK(BenchHasRepeatAutoVec<GenerateLowEntropySequence>)->Apply(ApplyHasRepeatArgs);
BENCHMARK(BenchHasRepeatAutoVec<GenerateMediumEntropySequence>)->Apply(ApplyHasRepeatArgs);
BENCHMARK(BenchHasRepeatAutoVec<GenerateHighEntropySequence>)->Apply(ApplyHasRepeatArgs);

// NOLINTEND(cert-err58-cpp, cppcoreguidelines-owning-memory, readability-identifier-length, misc-use-anonymous-namespace)
