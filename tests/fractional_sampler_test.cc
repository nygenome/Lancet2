#include "lancet2/fractional_sampler.h"

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

TEST_CASE("Can randomly downsample targeted fractions", "[lancet2::FractionalSampler]") {
  SECTION("sampling 0.3 fraction within 5% approximation") {
    lancet2::FractionalSampler sampler1{0.3};
    std::size_t hits1 = 0;
    for (std::size_t idx = 0; idx < 1e4; ++idx) {
      if (sampler1.ShouldSample()) hits1 += 1;
    }
    CHECK_THAT(hits1, Catch::Matchers::WithinRel(3e3, 0.05));
  }

  SECTION("sampling 0.85 fraction within 3% approximation") {
    lancet2::FractionalSampler sampler2{0.85};
    std::size_t hits2 = 0;
    for (std::size_t idx = 0; idx < 1e6; ++idx) {
      if (sampler2.ShouldSample()) hits2 += 1;
    }
    CHECK_THAT(hits2, Catch::Matchers::WithinRel(8.5e5, 0.03));
  }

  SECTION("samples all data with 1.0 fraction") {
    lancet2::FractionalSampler sampler3{1.0};
    std::size_t hits3 = 0;
    for (std::size_t idx = 0; idx < 1e6; ++idx) {
      if (sampler3.ShouldSample()) hits3 += 1;
    }
    CHECK(hits3 == 1e6);
  }
}
