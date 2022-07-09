#include "lancet2/fasta_reader.h"

#include "absl/strings/str_format.h"
#include "catch2/catch_test_macros.hpp"
#include "generated/test_config.h"

TEST_CASE("Can randomly access data from a FASTA file", "[lancet2::hts::FastaReader]") {
  const std::string faPath = absl::StrFormat("%s/human_g1k_v37.1_1_90000000.fa.gz", TEST_DATA_DIR);
  lancet2::FastaReader rdr{faPath};

  const auto regionSeq = rdr.GetRegionSeq(lancet2::GenomicRegion{"1", static_cast<u32>(80e6), static_cast<u32>(85e6)});
  REQUIRE(regionSeq.ok());
  CHECK(regionSeq.value().length() == 5000001);

  const auto ctgSeq = rdr.GetContigSeq("1");
  REQUIRE(ctgSeq.ok());
  CHECK(ctgSeq.value().length() == 90000000);

  const auto ctgs = rdr.GetContigs();
  REQUIRE(ctgs.size() == 1);
  CHECK(ctgs[0].contigName == "1");
  CHECK(ctgs[0].contigLen == 90000000);

  const auto ctgIdxMap = rdr.GetContigIndexMap();
  REQUIRE(ctgIdxMap.contains("1"));
  CHECK(ctgIdxMap.at("1") == 0);
  CHECK(ctgIdxMap.size() == 1);

  const auto ctgIdx = rdr.GetContigIndex("1");
  CHECK(ctgIdx.ok());
  CHECK(ctgIdx.value() == 0);

  const auto ctgLen = rdr.GetContigLength("1");
  CHECK(ctgLen.ok());
  CHECK(ctgLen.value() == 90000000);
}
