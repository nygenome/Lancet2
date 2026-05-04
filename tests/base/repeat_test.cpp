#include "lancet/base/repeat.h"

#include "lancet/base/types.h"

#include "catch_amalgamated.hpp"

#include <array>
#include <random>
#include <string>
#include <string_view>

namespace lancet::base::tests {

namespace {

/// Generate a 5000bp pseudo-random DNA sequence for Hamming distance fuzz
/// testing. Caller passes a const literal seed so the same TEST_CASE produces
/// the same sequence on every run — `std::random_device` is forbidden in tests
/// (per the project's determinism convention) because two failing runs would
/// be impossible to compare. Two distinct seeds in the same TEST_CASE produce
/// two distinct sequences (probabilistically) so the per-iteration comparison
/// `HammingDist(s, t) != 0` still exercises the "different inputs" path.
inline auto GenerateRandomDnaSequence(u64 seed) -> std::string {
  static constexpr std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};

  std::mt19937_64 generator(seed);

  static constexpr usize SEQ_LENGTH = 5000;
  std::uniform_int_distribution<usize> base_chooser(0, 3);
  std::string result;
  result.reserve(SEQ_LENGTH);

  for (usize idx = 0; idx < SEQ_LENGTH; ++idx) {
    result.push_back(BASES.at(base_chooser(generator)));
  }

  return result;
}

}  // namespace

// ============================================================================
//  Hamming distance: randomized correctness
// ============================================================================

// Catch2 SECTION fan-out inflates clang-tidy's cognitive-complexity metric beyond the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Can calculate hamming distance correctly for random strings",
          "[lancet][base][HammingDist]") {
  // Pinned base seed: every iteration derives a distinct seed by adding the
  // loop index, so the 1000 iterations produce 1000 deterministic sequences.
  // Re-running this TEST_CASE on a different machine (or after an unrelated
  // refactor) reproduces the same sequence pair on every iteration.
  static constexpr u64 BASE_SEED = 0x5E'ED'5E'ED'5E'ED'5E'EDULL;
  static constexpr usize NUM_ITERATIONS = 1000;

  for (usize idx = 0; idx <= NUM_ITERATIONS; ++idx) {
    auto const result = GenerateRandomDnaSequence(BASE_SEED + idx);
    // Offset the second seed far enough that the two streams diverge at the
    // very first base — keeps `HammingDist != 0` deterministically true.
    auto const other = GenerateRandomDnaSequence(BASE_SEED + idx + 0x10'00'00'00ULL);

    REQUIRE(HammingDist(result, result) == 0);
    REQUIRE(HammingDist(result, other) != 0);
  }
}

// ============================================================================
//  Hamming distance: small known-answer tests
// ============================================================================

TEST_CASE("Can calculate hamming distance correctly for small test",
          "[lancet][base][HammingDist]") {
  std::string_view const test = "aaaa";
  std::string_view const diff_a = "abaa";
  std::string_view const diff_b = "aaba";

  REQUIRE(HammingDist(test, test) == 0);
  REQUIRE(HammingDist(test, diff_a) == 1);
  REQUIRE(HammingDist(test, diff_b) == 1);
  REQUIRE(HammingDist(diff_a, diff_b) == 2);
}

}  // namespace lancet::base::tests
