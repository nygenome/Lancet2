#include <array>
#include <random>
#include <string>

#include "benchmark/benchmark.h"
#include "lancet/base/repeat.h"
#include "lancet/base/types.h"

[[nodiscard]] inline auto GenerateRandomDnaSequence(const usize seq_len) -> std::string {
  static constexpr std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};

  std::random_device device;
  std::mt19937_64 generator(device());

  std::uniform_int_distribution<usize> base_chooser(0, 3);
  std::string result(seq_len, 'N');

  for (usize idx = 0; idx < seq_len; ++idx) {
    result[idx] = BASES.at(base_chooser(generator));
  }

  return result;
}

namespace {

void BenchHammingNaive(benchmark::State& state) {
  const std::string first = GenerateRandomDnaSequence(static_cast<usize>(state.range(0)));
  const std::string second = GenerateRandomDnaSequence(static_cast<usize>(state.range(0)));

  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    auto result = HammingDistNaive(first, second);
    benchmark::DoNotOptimize(result);
  }
}

void BenchHamming64(benchmark::State& state) {
  const std::string first = GenerateRandomDnaSequence(static_cast<usize>(state.range(0)));
  const std::string second = GenerateRandomDnaSequence(static_cast<usize>(state.range(0)));

  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    auto result = HammingDistWord64(first, second);
    benchmark::DoNotOptimize(result);
  }
}

}  // namespace

// NOLINTBEGIN
BENCHMARK(BenchHammingNaive)->DenseRange(11, 101, 4);
BENCHMARK(BenchHamming64)->DenseRange(11, 101, 4);

BENCHMARK(BenchHammingNaive)->RangeMultiplier(2)->Range(2 << 2, 2 << 10);
BENCHMARK(BenchHamming64)->RangeMultiplier(2)->Range(2 << 2, 2 << 10);
// NOLINTEND
