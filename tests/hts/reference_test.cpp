#include "lancet/hts/reference.h"

#include <filesystem>
#include <optional>

#include "absl/status/status.h"
#include "catch_amalgamated.hpp"
#include "lancet/base/types.h"
#include "lancet_test_config.h"

using namespace lancet::hts;

TEST_CASE("Reference::Reference()", "[lancet][hts][Reference]") {
  const auto ref_path = MakePath(FULL_DATA_DIR, GRCH38_REF_NAME);
  REQUIRE(std::filesystem::exists(ref_path));
  REQUIRE_NOTHROW(Reference(ref_path));
  CHECK(Reference(ref_path).FastaPath() == ref_path);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reference::ListChroms()", "[lancet][hts][Reference]") {
  const Reference ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));
  const auto chromosomes = ref.ListChroms();
  static constexpr usize NUM_AUTOSOMES_XY = 24;

  SECTION("Reference must contain the expected number of chromosomes") { CHECK(chromosomes.size() == 3366); }

  SECTION("Reference must list the correct information for chromosomes") {
    for (usize idx = 0; idx < NUM_AUTOSOMES_XY; ++idx) {
      CHECK(chromosomes.at(idx).Name() == GRCH38_NAMES_AUTOSOMES_XY.at(idx));
      CHECK(chromosomes.at(idx).Length() == GRCH38_LENGTHS_AUTOSOMES_XY.at(idx));
      CHECK(chromosomes.at(idx).Index() == idx);
    }
  }
}

TEST_CASE("Reference::FindChromByName(std::string_view)", "[lancet][hts][Reference]") {
  const Reference ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));

  SECTION("Must return ok status and Chrom when finding a chromosome present in the reference") {
    const auto result = ref.FindChromByName("chrUn_JTFH01001570v1_decoy");
    REQUIRE(result.ok());
    CHECK(result->Name() == "chrUn_JTFH01001570v1_decoy");
    CHECK(result->Length() == 1244);
    CHECK(result->Index() == 2412);
  }

  SECTION("Must return NotFound status when finding a chromosome not present in the reference") {
    const auto result = ref.FindChromByName("missing-non-existing-chrom");
    REQUIRE_FALSE(result.ok());
    CHECK(result.status().code() == absl::StatusCode::kNotFound);
  }
}

TEST_CASE("Reference::MakeRegion(const char*)", "[lancet][hts][Reference]") {
  const Reference ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));

  SECTION("Can make Reference::Region from {HLA-B*15:01:01:02N}") {
    const auto region = ref.MakeRegion("{HLA-B*15:01:01:02N}");
    static constexpr auto expected_seq =
        "ATGCGGGTCACGGCGCCCCGAACCGTCCTCCTGCTGCTCTCGGGAGCCCTGGCCCTGACCGAGACCTGGGCCGGTGAGTGCGGGGTCGGCAGGGAAATGGCCTCTGTGGG"
        "GAGGAGCGAGGGGACCGCAGGCGGGGGCGCAGGACCCGGGGAGCCGCGCCGGGAGGAGGGTCGGGCCCCTCCTCGCCCCCAGGCTCCCACTCCATGAGGTATTTCTACAC"
        "CGCCATGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACCCAGTTCGTGAGGTTCGACAGCGACGCCGCGAGTCCGAGGATGGCGC"
        "CCCGGGCGCCATGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCGGGAGACACAGATCTCCAAGACCAACACACAGACTTACCGAGAGAGCCTGCGGAACCTGCGCGGC"
        "TACTACAACCAGAGCGAGGCCGGGTCTCACACCCTCCAGAGGATGTACGGCTGCGACGTGGGGCCGGACGGGCGCCTCCTCCGCGGGCATGACCAGTCCGCCTACGACGG"
        "CAAGGATTACATCGCCCTGAACGAGGACCTGAGCTCCTGGACCGCGGCGGACACGGCGGCTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTGAGGCGGAGCAGTGGA"
        "GAGCCTACCTGGAGGGCCTGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGGACCCCCCAAAGACACATGTGACCCACCACCCC"
        "ATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGCGAGGACCAAACTCAGGACACCGAGCT"
        "TGTGGAGACCAGACCAGCAGGAGATAGAACCTTCCAGAAGTGGGCAGCTGTGGTGGTGCCTTCTGGAGAAGAGCAGAGATACACATGCCATGTACAGCATGAGGGGCTGC"
        "CGAAGCCCCTCACCCTGAGATGGGAGCCATCTTCCCAGTCCACCATCCCCATCGTGGGCATTGTTGCTGGCCTGGCTGTCCTAGCAGTTGTGGTCATCGGAGCTGTGGTC"
        "GCTACTGTGATGTGTAGGAGGAAGAGCTCAGGTGGAAAAGGAGGGAGCTACTCTCAGGCTGCGTCCAGCGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGA";

    CHECK(region.ChromName() == "HLA-B*15:01:01:02N");
    CHECK(region.ChromIndex() == 3003);
    CHECK(region.StartPos1() == 1);
    CHECK(region.EndPos1() == region.Length());
    CHECK(region.Length() == 1208);
    CHECK(region.SeqView() == expected_seq);
  }

  SECTION("Can make Reference::Region from {HLA-A*01:01:01:01}:100-110") {
    const auto region = ref.MakeRegion("{HLA-A*01:01:01:01}:100-110");
    static constexpr auto expected_seq = "GGGGATTCCCC";
    CHECK(region.ChromName() == "HLA-A*01:01:01:01");
    CHECK(region.ChromIndex() == 2841);
    CHECK(region.StartPos1() == 100);
    CHECK(region.EndPos1() == 110);
    CHECK(region.Length() == 11);
    CHECK(region.SeqView() == expected_seq);
  }

  SECTION("Can make Reference::Region from {HLA-A*01:01:01:01}:3434-") {
    const auto region = ref.MakeRegion("{HLA-A*01:01:01:01}:3434-");
    static constexpr auto expected_seq = "CTGAGGTGTCTCCATCTCTGTCTCAACTTCATGGTGCACTGAGCTGTAACTTCTTCCTTCCCTATTAAAA";
    CHECK(region.ChromName() == "HLA-A*01:01:01:01");
    CHECK(region.ChromIndex() == 2841);
    CHECK(region.StartPos1() == 3434);
    CHECK(region.EndPos1() == 3503);
    CHECK(region.Length() == 70);
    CHECK(region.SeqView() == expected_seq);
  }

  SECTION("Can make Reference::Region from {HLA-A*01:01:01:01}:-5") {
    const auto region = ref.MakeRegion("{HLA-A*01:01:01:01}:-5");
    static constexpr auto expected_seq = "CAGGA";
    CHECK(region.ChromName() == "HLA-A*01:01:01:01");
    CHECK(region.ChromIndex() == 2841);
    CHECK(region.StartPos1() == 1);
    CHECK(region.EndPos1() == 5);
    CHECK(region.Length() == 5);
    CHECK(region.SeqView() == expected_seq);
  }

  SECTION("Can make Reference::Region chr18:48343-48343") {
    const auto region = ref.MakeRegion("chr18:48343-48343");
    static constexpr auto expected_seq = "C";
    CHECK(region.ChromName() == "chr18");
    CHECK(region.ChromIndex() == 17);
    CHECK(region.StartPos1() == 48343);
    CHECK(region.EndPos1() == 48343);
    CHECK(region.Length() == 1);
    CHECK(region.SeqView() == expected_seq);
  }
}

TEST_CASE("Reference::MakeRegion(const std::string&, const OneBasedClosedOptional&)", "[lancet][hts][Reference]") {
  const Reference ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));

  SECTION("Can make Reference::Region from (\"HLA-B*15:01:01:02N\", {})") {
    const auto region = ref.MakeRegion("HLA-B*15:01:01:02N", {});
    static constexpr auto expected_seq =
        "ATGCGGGTCACGGCGCCCCGAACCGTCCTCCTGCTGCTCTCGGGAGCCCTGGCCCTGACCGAGACCTGGGCCGGTGAGTGCGGGGTCGGCAGGGAAATGGCCTCTGTGGG"
        "GAGGAGCGAGGGGACCGCAGGCGGGGGCGCAGGACCCGGGGAGCCGCGCCGGGAGGAGGGTCGGGCCCCTCCTCGCCCCCAGGCTCCCACTCCATGAGGTATTTCTACAC"
        "CGCCATGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACCCAGTTCGTGAGGTTCGACAGCGACGCCGCGAGTCCGAGGATGGCGC"
        "CCCGGGCGCCATGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCGGGAGACACAGATCTCCAAGACCAACACACAGACTTACCGAGAGAGCCTGCGGAACCTGCGCGGC"
        "TACTACAACCAGAGCGAGGCCGGGTCTCACACCCTCCAGAGGATGTACGGCTGCGACGTGGGGCCGGACGGGCGCCTCCTCCGCGGGCATGACCAGTCCGCCTACGACGG"
        "CAAGGATTACATCGCCCTGAACGAGGACCTGAGCTCCTGGACCGCGGCGGACACGGCGGCTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTGAGGCGGAGCAGTGGA"
        "GAGCCTACCTGGAGGGCCTGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGGACCCCCCAAAGACACATGTGACCCACCACCCC"
        "ATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGCGAGGACCAAACTCAGGACACCGAGCT"
        "TGTGGAGACCAGACCAGCAGGAGATAGAACCTTCCAGAAGTGGGCAGCTGTGGTGGTGCCTTCTGGAGAAGAGCAGAGATACACATGCCATGTACAGCATGAGGGGCTGC"
        "CGAAGCCCCTCACCCTGAGATGGGAGCCATCTTCCCAGTCCACCATCCCCATCGTGGGCATTGTTGCTGGCCTGGCTGTCCTAGCAGTTGTGGTCATCGGAGCTGTGGTC"
        "GCTACTGTGATGTGTAGGAGGAAGAGCTCAGGTGGAAAAGGAGGGAGCTACTCTCAGGCTGCGTCCAGCGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGA";

    CHECK(region.ChromName() == "HLA-B*15:01:01:02N");
    CHECK(region.ChromIndex() == 3003);
    CHECK(region.StartPos1() == 1);
    CHECK(region.EndPos1() == region.Length());
    CHECK(region.Length() == 1208);
    CHECK(region.SeqView() == expected_seq);
  }

  SECTION("Can make Reference::Region from (\"HLA-A*01:01:01:01\", {100, 110})") {
    const auto region = ref.MakeRegion("HLA-A*01:01:01:01", {100, 110});
    static constexpr auto expected_seq = "GGGGATTCCCC";
    CHECK(region.ChromName() == "HLA-A*01:01:01:01");
    CHECK(region.ChromIndex() == 2841);
    CHECK(region.StartPos1() == 100);
    CHECK(region.EndPos1() == 110);
    CHECK(region.Length() == 11);
    CHECK(region.SeqView() == expected_seq);
  }

  SECTION("Can make Reference::Region from (\"HLA-A*01:01:01:01\", {3434, std::nullopt})") {
    const auto region = ref.MakeRegion("HLA-A*01:01:01:01", {3434, std::nullopt});
    static constexpr auto expected_seq = "CTGAGGTGTCTCCATCTCTGTCTCAACTTCATGGTGCACTGAGCTGTAACTTCTTCCTTCCCTATTAAAA";
    CHECK(region.ChromName() == "HLA-A*01:01:01:01");
    CHECK(region.ChromIndex() == 2841);
    CHECK(region.StartPos1() == 3434);
    CHECK(region.EndPos1() == 3503);
    CHECK(region.Length() == 70);
    CHECK(region.SeqView() == expected_seq);
  }

  SECTION("Can make Reference::Region from (\"HLA-A*01:01:01:01\", {std::nullopt, 5})") {
    const auto region = ref.MakeRegion("HLA-A*01:01:01:01", {std::nullopt, 5});
    static constexpr auto expected_seq = "CAGGA";
    CHECK(region.ChromName() == "HLA-A*01:01:01:01");
    CHECK(region.ChromIndex() == 2841);
    CHECK(region.StartPos1() == 1);
    CHECK(region.EndPos1() == 5);
    CHECK(region.Length() == 5);
    CHECK(region.SeqView() == expected_seq);
  }

  SECTION("Can make Reference::Region from (\"chr18\", {48343, 48343})") {
    const auto region = ref.MakeRegion("chr18", {48343, 48343});
    static constexpr auto expected_seq = "C";
    CHECK(region.ChromName() == "chr18");
    CHECK(region.ChromIndex() == 17);
    CHECK(region.StartPos1() == 48343);
    CHECK(region.EndPos1() == 48343);
    CHECK(region.Length() == 1);
    CHECK(region.SeqView() == expected_seq);
  }
}
