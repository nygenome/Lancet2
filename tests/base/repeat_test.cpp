#include "lancet/base/repeat.h"

#include "lancet/base/types.h"

#include "catch_amalgamated.hpp"

#include <array>
#include <random>
#include <string>
#include <string_view>

using lancet::base::HammingDist;

namespace {

/// Generate a 5000bp random DNA sequence for Hamming distance fuzz testing.
inline auto GenerateRandomDnaSequence() -> std::string {
  static constexpr std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};

  std::random_device device;
  std::mt19937_64 generator(device());

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
  static constexpr usize NUM_ITERATIONS = 1000;

  for (usize idx = 0; idx <= NUM_ITERATIONS; ++idx) {
    auto const result = GenerateRandomDnaSequence();
    auto const other = GenerateRandomDnaSequence();

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
