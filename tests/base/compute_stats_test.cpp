#include "lancet/base/compute_stats.h"

#include "lancet/base/types.h"

#include "absl/random/distributions.h"
#include "absl/types/span.h"
#include "catch_amalgamated.hpp"

#include <array>
#include <random>
#include <vector>

#include <cmath>

namespace lancet::base::tests {

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  OnlineStats — Welford accumulator                                       ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("OnlineStats reports zero mean and variance when empty", "[lancet][base][OnlineStats]") {
  OnlineStats const stats;
  CHECK(stats.IsEmpty());
  CHECK(stats.Count() == 0);
  CHECK(stats.Mean() == Catch::Approx(0.0).margin(1e-12));
  CHECK(stats.Variance() == Catch::Approx(0.0).margin(1e-12));
  CHECK(stats.StdDev() == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE(
    "OnlineStats Add tracks Mean / Variance / StdDev on the canonical {2,4,4,4,5,5,7,9} sample",
    "[lancet][base][OnlineStats]") {
  // Hand-checked reference: μ = 5, σ² (Bessel-corrected, n=8, divisor 7) = 32/7,
  // σ = √(32/7). The {2,4,4,4,5,5,7,9} sample is the standard Welford textbook
  // example: small enough to verify by hand, large enough to exercise both the
  // running-mean update and the m2 accumulator.
  OnlineStats stats;
  for (auto const value : std::array<i32, 8>{2, 4, 4, 4, 5, 5, 7, 9}) {
    stats.Add(value);
  }
  CHECK_FALSE(stats.IsEmpty());
  CHECK(stats.Count() == 8);
  CHECK(stats.Mean() == Catch::Approx(5.0).margin(1e-12));
  CHECK(stats.Variance() == Catch::Approx(32.0 / 7.0).margin(1e-12));
  CHECK(stats.StdDev() == Catch::Approx(std::sqrt(32.0 / 7.0)).margin(1e-12));
}

TEST_CASE("OnlineStats Variance is 0 when only one value has been added",
          "[lancet][base][OnlineStats]") {
  // Bessel's correction divides by (n − 1); with n = 1 the implementation
  // returns 0 to avoid division by zero. Mean of a one-element stream is
  // the value itself.
  OnlineStats stats;
  stats.Add(42.0);
  CHECK(stats.Count() == 1);
  CHECK(stats.Mean() == Catch::Approx(42.0).margin(1e-12));
  CHECK(stats.Variance() == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("OnlineStats Clear returns the accumulator to its initial state",
          "[lancet][base][OnlineStats]") {
  OnlineStats stats;
  stats.Add(10.0);
  stats.Add(20.0);
  stats.Clear();
  CHECK(stats.IsEmpty());
  CHECK(stats.Count() == 0);
  CHECK(stats.Mean() == Catch::Approx(0.0).margin(1e-12));
  CHECK(stats.Variance() == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("OnlineStats Merge of three accumulators matches a single-pass run",
          "[lancet][base][OnlineStats]") {
  // Per-thread accumulators must merge into the same statistics a single
  // sequential pass would have produced. This validates Chan et al.'s
  // parallel m2 update m2_combined = m2_1 + m2_2 + δ²·n_1·n_2/(n_1+n_2).
  std::array<f64, 9> const values{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

  OnlineStats sequential;
  for (auto const value : values) sequential.Add(value);

  // Three shards of unequal size to exercise the merge formula in both
  // directions (small into large, large into small).
  OnlineStats shard_a;  // values 0..1
  OnlineStats shard_b;  // values 2..5
  OnlineStats shard_c;  // values 6..8
  for (usize idx = 0; idx < 2; ++idx) shard_a.Add(values[idx]);
  for (usize idx = 2; idx < 6; ++idx) shard_b.Add(values[idx]);
  for (usize idx = 6; idx < 9; ++idx) shard_c.Add(values[idx]);

  OnlineStats merged;
  merged.Merge(shard_a);
  merged.Merge(shard_b);
  merged.Merge(shard_c);

  CHECK(merged.Count() == sequential.Count());
  CHECK(merged.Mean() == Catch::Approx(sequential.Mean()).margin(1e-12));
  CHECK(merged.Variance() == Catch::Approx(sequential.Variance()).margin(1e-12));
  CHECK(merged == sequential);
}

TEST_CASE("OnlineStats Merge with an empty accumulator is a no-op", "[lancet][base][OnlineStats]") {
  OnlineStats populated;
  for (auto const value : std::array<f64, 4>{1.0, 2.0, 3.0, 4.0}) populated.Add(value);

  OnlineStats merged = populated;
  OnlineStats const empty_other;
  merged.Merge(empty_other);

  CHECK(merged == populated);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  Free Mean / Median / Minimum on integer spans                           ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("Mean of an empty span is 0.0", "[lancet][base][Mean]") {
  std::vector<i32> const data;
  CHECK(Mean(absl::MakeConstSpan(data)) == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("Mean of a single-element span returns the element", "[lancet][base][Mean]") {
  std::vector<i32> const data{7};
  CHECK(Mean(absl::MakeConstSpan(data)) == Catch::Approx(7.0).margin(1e-12));
}

TEST_CASE("Mean of {1..10} is 5.5", "[lancet][base][Mean]") {
  std::vector<i32> const data{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  CHECK(Mean(absl::MakeConstSpan(data)) == Catch::Approx(5.5).margin(1e-12));
}

TEST_CASE("Median of an empty span is 0", "[lancet][base][Median]") {
  std::vector<i32> const data;
  CHECK(Median(absl::MakeConstSpan(data)) == 0);
}

TEST_CASE("Median of a single-element span returns the element", "[lancet][base][Median]") {
  std::vector<i32> const data{42};
  CHECK(Median(absl::MakeConstSpan(data)) == 42);
}

TEST_CASE("Median of an odd-length sorted span returns the middle element",
          "[lancet][base][Median]") {
  std::vector<i32> const data{1, 3, 5, 7, 9};
  CHECK(Median(absl::MakeConstSpan(data)) == 5);
}

TEST_CASE("Median of an even-length sorted span returns the average of the two middle elements",
          "[lancet][base][Median]") {
  std::vector<i32> const data{1, 2, 3, 4};
  // Truncating integer division: (2 + 3) / 2 == 2.
  CHECK(Median(absl::MakeConstSpan(data)) == 2);
}

TEST_CASE("Median is order-independent (unsorted input)", "[lancet][base][Median]") {
  std::vector<i32> const data{9, 1, 5, 3, 7};
  CHECK(Median(absl::MakeConstSpan(data)) == 5);
}

TEST_CASE("Minimum of an empty span is 0", "[lancet][base][Minimum]") {
  std::vector<i32> const data;
  CHECK(Minimum(absl::MakeConstSpan(data)) == 0);
}

TEST_CASE("Minimum returns the smallest element on a non-empty span", "[lancet][base][Minimum]") {
  std::vector<i32> const data{4, 2, 9, 1, 7};
  CHECK(Minimum(absl::MakeConstSpan(data)) == 1);
}

TEST_CASE("Minimum tolerates negative values", "[lancet][base][Minimum]") {
  std::vector<i32> const data{4, -5, 2, -10, 7};
  CHECK(Minimum(absl::MakeConstSpan(data)) == -10);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  Property: OnlineStats.Variance() matches the two-pass formula           ║
// ║                                                                          ║
// ║  Welford's recurrence (m1, m2) is supposed to compute the same variance  ║
// ║  as the textbook two-pass formula                                        ║
// ║      Var = (1/(n-1)) * Σ(xᵢ − μ)²                                        ║
// ║  but with better numerical stability when the mean is large relative to  ║
// ║  the spread. Both formulas should agree on samples where the naive       ║
// ║  computation doesn't lose precision — that's exactly the input class     ║
// ║  this test exercises (uniform on [0, 100]).                              ║
// ╚══════════════════════════════════════════════════════════════════════════╝

// Catch2's per-iteration generation drives the cognitive-complexity metric over the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("OnlineStats.Variance matches the two-pass formula on randomly generated samples",
          "[lancet][base][OnlineStats]") {
  // Pinned seed: deterministic across platforms and reruns. A regression
  // that broke Welford's recurrence in some specific input subset would
  // surface here on the same iteration every run.
  static constexpr u64 BASE_SEED = 0x5E'ED'5E'ED'5E'ED'5E'EDULL;
  static constexpr usize NUM_PROPERTY_ITERATIONS = 200;
  static constexpr usize SAMPLE_SIZE = 50;

  // Const-literal seed is the project's documented determinism convention
  // (see test_style.md / §A.9). The clang-tidy check is conservatively
  // designed for production code; in tests, predictability is exactly
  // what we want.
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937_64 generator(BASE_SEED);

  for (usize iter = 0; iter < NUM_PROPERTY_ITERATIONS; ++iter) {
    std::vector<f64> samples(SAMPLE_SIZE);
    for (auto& sample : samples) {
      sample = absl::Uniform<f64>(absl::IntervalClosed, generator, 0.0, 100.0);
    }

    // Welford via the in-tree implementation.
    OnlineStats welford;
    for (auto const sample : samples) welford.Add(sample);

    // Two-pass reference: first pass for mean, second pass for variance.
    f64 sum = 0.0;
    for (auto const sample : samples) sum += sample;
    f64 const mean = sum / static_cast<f64>(samples.size());

    f64 sum_sq_dev = 0.0;
    for (auto const sample : samples) {
      f64 const dev = sample - mean;
      sum_sq_dev += dev * dev;
    }
    f64 const two_pass_variance = sum_sq_dev / static_cast<f64>(samples.size() - 1);

    INFO("iter=" << iter);
    // 1e-9 absorbs the f64 rounding accumulated through Welford's update
    // path vs the two summations in the reference. Anything tighter risks
    // false-failing on platform-specific FP rounding.
    CHECK(welford.Variance() == Catch::Approx(two_pass_variance).margin(1e-9));
  }
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  Property: OnlineStats.Merge is associative                              ║
// ║                                                                          ║
// ║  Chan et al.'s parallel-merge formula must produce the same combined     ║
// ║  accumulator regardless of the order in which three shards are merged:   ║
// ║      ((A ∪ B) ∪ C) == (A ∪ (B ∪ C)) == ((A ∪ C) ∪ B)                     ║
// ║  Associativity is the foundation that lets the runtime merge per-thread  ║
// ║  shards in any tree shape without changing the result.                   ║
// ╚══════════════════════════════════════════════════════════════════════════╝

// Catch2's per-iteration generation drives the cognitive-complexity metric over the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("OnlineStats.Merge is associative across three accumulators",
          "[lancet][base][OnlineStats]") {
  // Pinned seed for deterministic input across runs.
  static constexpr u64 BASE_SEED = 0x5E'ED'5E'ED'5E'ED'5E'EDULL;
  static constexpr usize NUM_PROPERTY_ITERATIONS = 100;

  // Const-literal seed is the project's documented determinism convention
  // (see test_style.md / §A.9). The clang-tidy check is conservatively
  // designed for production code; in tests, predictability is exactly
  // what we want.
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937_64 generator(BASE_SEED);

  for (usize iter = 0; iter < NUM_PROPERTY_ITERATIONS; ++iter) {
    auto const fill_shard = [&](OnlineStats& shard) {
      auto const shard_size = absl::Uniform<usize>(absl::IntervalClosed, generator, 1, 30);
      for (usize sample_idx = 0; sample_idx < shard_size; ++sample_idx) {
        shard.Add(absl::Uniform<f64>(absl::IntervalClosed, generator, 0.0, 100.0));
      }
    };

    OnlineStats shard_a;
    OnlineStats shard_b;
    OnlineStats shard_c;
    fill_shard(shard_a);
    fill_shard(shard_b);
    fill_shard(shard_c);

    // Three associativity parens: ((A∪B)∪C), (A∪(B∪C)), ((A∪C)∪B).
    OnlineStats merged_left = shard_a;
    merged_left.Merge(shard_b);
    merged_left.Merge(shard_c);

    OnlineStats merged_right_inner = shard_b;
    merged_right_inner.Merge(shard_c);
    OnlineStats merged_right = shard_a;
    merged_right.Merge(merged_right_inner);

    OnlineStats merged_swap = shard_a;
    merged_swap.Merge(shard_c);
    merged_swap.Merge(shard_b);

    // Associativity is exact in real arithmetic but not bit-exact in f64
    // because Chan's merge formula accumulates rounding differently when
    // merge order shifts the relative magnitudes of the (n_1, n_2)
    // weights. `OnlineStats::operator==` uses `epsilon ≈ 2.22e-16` which
    // is too tight for these accumulated errors; we compare Mean and
    // Variance directly with a looser tolerance that's still well below
    // any biologically-meaningful difference.
    constexpr f64 ASSOCIATIVITY_TOLERANCE = 1e-9;
    INFO("iter=" << iter);
    CHECK(merged_left.Count() == merged_right.Count());
    CHECK(merged_left.Count() == merged_swap.Count());
    CHECK(merged_left.Mean() == Catch::Approx(merged_right.Mean()).margin(ASSOCIATIVITY_TOLERANCE));
    CHECK(merged_left.Mean() == Catch::Approx(merged_swap.Mean()).margin(ASSOCIATIVITY_TOLERANCE));
    CHECK(merged_left.Variance() ==
          Catch::Approx(merged_right.Variance()).margin(ASSOCIATIVITY_TOLERANCE));
    CHECK(merged_left.Variance() ==
          Catch::Approx(merged_swap.Variance()).margin(ASSOCIATIVITY_TOLERANCE));
  }
}

}  // namespace lancet::base::tests
