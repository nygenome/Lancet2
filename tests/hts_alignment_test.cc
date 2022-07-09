#include "lancet2/hts_alignment.h"

#include "absl/strings/str_format.h"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_vector.hpp"
#include "generated/test_config.h"
#include "lancet2/assert_macro.h"
#include "lancet2/hts_reader.h"

static inline auto GetFirstAlignment() -> lancet2::HtsAlignment {
  const auto bamPath = absl::StrFormat("%s/tumor.bam", TEST_DATA_DIR);
  const auto refPath = absl::StrFormat("%s/human_g1k_v37.1_1_90000000.fa.gz", TEST_DATA_DIR);

  lancet2::HtsReader rdr(bamPath, refPath);
  static lancet2::HtsAlignment aln;
  const auto state = rdr.GetNextAlignment(&aln, {"XT", "XA", "AS", "XS"});
  LANCET_ASSERT(state == lancet2::HtsReader::IteratorState::VALID);
  return aln;
}

static inline auto MakeBQVector(const std::string& bq) -> std::vector<int> {
  std::vector<int> result;
  result.reserve(bq.size());
  for (const auto& c : bq) {
    result.emplace_back(static_cast<int>(c));
  }

  return result;
}

TEST_CASE("Can fetch and query data from an alignment", "[lancet2::HtsAlignment]") {
  const auto al = GetFirstAlignment();
  REQUIRE(al.GetLength() == 151);

  CHECK(al.GetReadSequence() ==
        "GCTGCAATATGATTTTTTAAAGTCATTTCTTGATCATTAAACACTTTAATTTGCTGCATTTCAAATCATACTATCTGTCCCAGGTTAGATGTAGTTCATGCTGAGAAATG"
        "TGGGTAAAGAATTCAAGTCCATTTTAAAATTCTCCAATTCT");

  CHECK_THAT(MakeBQVector(al.GetReadQuality()),
             Catch::Matchers::Equals(std::vector<int>{
                 26, 24, 25, 26, 23, 29, 30, 31, 31, 31, 29, 30, 31, 31, 31, 30, 30, 30, 29, 29, 31, 29, 30, 30, 25, 30,
                 29, 30, 24, 26, 30, 29, 30, 30, 30, 30, 30, 30, 28, 29, 25, 30, 30, 30, 30, 30, 30, 29, 30, 27, 27, 30,
                 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 29, 30,
                 30, 30, 30, 31, 30, 29, 30, 30, 31, 29, 30, 30, 29, 30, 31, 29, 30, 30, 30, 30, 30, 30, 30, 30, 29, 31,
                 29, 29, 29, 30, 30, 29, 30, 30, 30, 29, 30, 29, 29, 31, 29, 29, 30, 30, 29, 29, 28, 30, 28, 29, 29, 29,
                 29, 29, 29, 29, 29, 28, 28, 28, 30, 29, 29, 29, 29, 29, 30, 29, 30, 30, 28, 28, 30}));

  CHECK(al.GetReadName() == "E00217:32:H25FNCCXX:1:2117:32170:31758");
  CHECK(al.GetContigName() == "1");
  CHECK(al.GetMateContigName() == "1");
  CHECK(al.GetStartPos0() == 82959850);
  CHECK(al.GetEndPos0() == 82960001);
  CHECK(al.GetMateStartPos0() == 82959683);
  CHECK(al.GetMappingQual() == 60);

  CHECK(al.IsPaired());
  CHECK(al.IsProperPair());
  CHECK(al.IsReverseStrand());
  CHECK(al.IsRead1());
  CHECK(al.IsPrimary());

  CHECK_FALSE(al.IsDuplicate());
  CHECK_FALSE(al.IsSupplementary());
  CHECK_FALSE(al.IsSecondary());
  CHECK_FALSE(al.IsQcFailed());
  CHECK_FALSE(al.IsUnmapped());
  CHECK_FALSE(al.IsMateUnmapped());
  CHECK_FALSE(al.IsMateReverseStrand());
  CHECK_FALSE(al.IsRead2());

  CHECK(lancet2::ToString(al.GetCigarData()) == "151M");
  CHECK_NOTHROW(al.GetMateRegion());
  CHECK(al.HasTag("AS"));
  CHECK(al.HasTag("XS"));
  CHECK_NOTHROW(al.GetTagData("AS"));

  CHECK(al.IsWithinRegion(lancet2::GenomicRegion{"1", 82000000, 89000000}));
  CHECK_FALSE(al.IsWithinRegion(lancet2::GenomicRegion{"10", 100, 1000}));

  std::vector<std::uint32_t> vals;
  CHECK_FALSE(al.GetSoftClips(&vals, &vals, &vals));

  SECTION("Can strip bases from both 5' and 3' ends to build ReadInfo") {
    // After stripping first 7 and last 33 bases, we end up with 111 bases
    const auto ri = al.BuildReadInfo(lancet2::SampleLabel::TUMOR, 31, 101);
    REQUIRE_FALSE(ri.IsEmpty());
    CHECK(ri.Length() == 111);
    CHECK_FALSE(ri.IsEmpty());

    CHECK(ri.sequence ==
          "TATGATTTTTTAAAGTCATTTCTTGATCATTAAACACTTTAATTTGCTGCATTTCAAATCATACTATCTGTCCCAGGTTAGATGTAGTTCATGCTGAGAAATG"
          "TGGGTAAA");

    CHECK_THAT(MakeBQVector(ri.quality),
               Catch::Matchers::Equals(std::vector<int>{
                   31, 31, 31, 29, 30, 31, 31, 31, 30, 30, 30, 29, 29, 31, 29, 30, 30, 25, 30, 29, 30, 24, 26,
                   30, 29, 30, 30, 30, 30, 30, 30, 28, 29, 25, 30, 30, 30, 30, 30, 30, 29, 30, 27, 27, 30, 30,
                   30, 30, 30, 30, 30, 30, 30, 30, 30, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
                   29, 30, 30, 30, 30, 31, 30, 29, 30, 30, 31, 29, 30, 30, 29, 30, 31, 29, 30, 30, 30, 30, 30,
                   30, 30, 30, 29, 31, 29, 29, 29, 30, 30, 29, 30, 30, 30, 29, 30, 29, 29, 31}));
  }

  SECTION("Can read min final len and skip building ReadInfo") {
    const auto result = al.BuildReadInfo(lancet2::SampleLabel::TUMOR, 31, 121);
    REQUIRE(result.IsEmpty());
  }

  SECTION("Can read min trim BQ and skip building ReadInfo") {
    const auto result = al.BuildReadInfo(lancet2::SampleLabel::TUMOR, 45, 100);
    REQUIRE(result.IsEmpty());
  }
}
