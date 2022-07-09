#include "lancet2/genomic_region.h"

#include "catch2/catch_test_macros.hpp"

TEST_CASE("Can build and get information from a hts region", "[lancet2::hts::HtsRegion]") {
  SECTION("hts region with chrom, start and end set") {
    lancet2::GenomicRegion r1{"chr1", 10, 100};
    CHECK(r1.GetChromName() == "chr1");
    CHECK(r1.GetStartPos1() == 10);
    CHECK(r1.GetEndPos1() == 100);
    CHECK(r1.GetLength() == 91);
    CHECK(r1.ToSamtoolsRegion() == "chr1:10-100");
  }

  SECTION("hts region with only chrom set") {
    lancet2::GenomicRegion r2{"chr10"};
    CHECK(r2.GetChromName() == "chr10");
    CHECK(r2.GetStartPos1() == -1);
    CHECK(r2.GetEndPos1() == -1);
    CHECK(r2.GetLength() == 0);
    CHECK(r2.ToSamtoolsRegion() == "chr10");
  }
}
