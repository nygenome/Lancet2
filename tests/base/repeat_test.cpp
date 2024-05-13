#include "lancet/base/repeat.h"

#include <array>
#include <random>
#include <string>
#include <string_view>

#include "catch_amalgamated.hpp"
#include "lancet/base/types.h"

namespace {

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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Can calculate hamming distance correctly for random strings", "[lancet][base][repeat]") {
  static constexpr usize NUM_ITERATIONS = 1000;

  for (usize idx = 0; idx <= NUM_ITERATIONS; ++idx) {
    const auto result = GenerateRandomDnaSequence();
    const auto other = GenerateRandomDnaSequence();

    SECTION("Naive method calculates correct distances") {
      REQUIRE(HammingDistNaive(result, result) == 0);
      REQUIRE(HammingDistNaive(result, other) != 0);
    }

    SECTION("64-bit word method calculates correct distances") {
      REQUIRE(HammingDistWord64(result, result) == 0);
      REQUIRE(HammingDistWord64(result, other) != 0);
    }
  }
}

TEST_CASE("Can calculate hamming distance correctly for small test", "[lancet][base][repeat]") {
  const std::string_view test = "aaaa";
  const std::string_view diff_a = "abaa";
  const std::string_view diff_b = "aaba";

  SECTION("Naive method calculates correct distances") {
    REQUIRE(HammingDistNaive(test, test) == 0);
    REQUIRE(HammingDistNaive(test, diff_a) == 1);
    REQUIRE(HammingDistNaive(test, diff_b) == 1);
    REQUIRE(HammingDistNaive(diff_a, diff_b) == 2);
  }

  SECTION("64-bit word method calculates correct distances") {
    REQUIRE(HammingDistWord64(test, test) == 0);
    REQUIRE(HammingDistWord64(test, diff_a) == 1);
    REQUIRE(HammingDistWord64(test, diff_b) == 1);
    REQUIRE(HammingDistWord64(diff_a, diff_b) == 2);
  }
}
