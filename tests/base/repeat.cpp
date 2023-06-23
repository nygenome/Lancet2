#include "lancet/base/repeat.h"

#include <string_view>

#include "catch_amalgamated.hpp"

TEST_CASE("Can calculcate hamming distance correctly", "[lancet][base][repeat]") {
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
