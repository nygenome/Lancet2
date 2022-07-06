#include "lancet2/genomic_region.h"

#include "catch2/catch_test_macros.hpp"

TEST_CASE("Can build and get information from a hts region", "[lancet2::hts::HtsRegion]") {
  SECTION("hts region with chrom, start and end set") {
    lancet2::GenomicRegion r1{"chr1", 10, 100};
    CHECK(r1.Chromosome() == "chr1");
    CHECK(r1.StartPosition1() == 10);
    CHECK(r1.EndPosition1() == 100);
    CHECK(r1.Length() == 91);
    CHECK(r1.ToRegionString() == "chr1:10-100");
  }

  SECTION("hts region with only chrom set") {
    lancet2::GenomicRegion r2{"chr10"};
    CHECK(r2.Chromosome() == "chr10");
    CHECK(r2.StartPosition1() == -1);
    CHECK(r2.EndPosition1() == -1);
    CHECK(r2.Length() == 0);
    CHECK(r2.ToRegionString() == "chr10");
  }
}
