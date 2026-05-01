#include "lancet/hts/extractor.h"

#include "lancet/hts/alignment.h"
#include "lancet/hts/reference.h"

#include "absl/types/span.h"
#include "catch_amalgamated.hpp"
#include "lancet_test_config.h"

#include <filesystem>
#include <initializer_list>
#include <iterator>
#include <stdexcept>

using lancet::hts::Alignment;
using lancet::hts::Extractor;
using lancet::hts::Reference;

TEST_CASE("Extractor::Extractor()", "[lancet][hts][Extractor]") {
  Reference const ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));
  auto const case_cram_path = MakePath(FULL_DATA_DIR, CASE_CRAM_NAME);
  auto const case_bam_path = MakePath(FULL_DATA_DIR, CASE_BAM_NAME);
  auto const ctrl_cram_path = MakePath(FULL_DATA_DIR, CTRL_CRAM_NAME);
  auto const ctrl_bam_path = MakePath(FULL_DATA_DIR, CTRL_BAM_NAME);

  REQUIRE(std::filesystem::exists(case_cram_path));
  REQUIRE(std::filesystem::exists(case_bam_path));
  REQUIRE(std::filesystem::exists(ctrl_cram_path));
  REQUIRE(std::filesystem::exists(ctrl_bam_path));

  CHECK_NOTHROW(Extractor(case_cram_path, ref, Alignment::Fields::CORE_QNAME));
  CHECK_NOTHROW(Extractor(case_bam_path, ref, Alignment::Fields::CORE_QNAME));
  CHECK_NOTHROW(Extractor(ctrl_cram_path, ref, Alignment::Fields::CORE_QNAME));
  CHECK_NOTHROW(Extractor(ctrl_bam_path, ref, Alignment::Fields::CORE_QNAME));
}

TEST_CASE("Extractor::SampleName()", "[lancet][hts][Extractor]") {
  Reference const ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));
  auto const case_cram_path = MakePath(FULL_DATA_DIR, CASE_CRAM_NAME);
  auto const case_bam_path = MakePath(FULL_DATA_DIR, CASE_BAM_NAME);

  Extractor const cram_extractor(case_cram_path, ref, Alignment::Fields::CORE_QNAME);
  Extractor const bam_extractor(case_bam_path, ref, Alignment::Fields::CORE_QNAME);

  CHECK(cram_extractor.SampleName() == bam_extractor.SampleName());
  CHECK(cram_extractor.SampleName() == "SRR7890893");
}

TEST_CASE("Extractor::ChromName(i32)", "[lancet][hts][Extractor]") {
  Reference const ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));
  auto const case_cram_path = MakePath(FULL_DATA_DIR, CASE_CRAM_NAME);
  auto const case_bam_path = MakePath(FULL_DATA_DIR, CASE_BAM_NAME);

  Extractor const cram_extractor(case_cram_path, ref, Alignment::Fields::CORE_QNAME);
  Extractor const bam_extractor(case_bam_path, ref, Alignment::Fields::CORE_QNAME);

  CHECK(cram_extractor.ChromName(10) == "chr11");
  CHECK(cram_extractor.ChromName(2841) == "HLA-A*01:01:01:01");
  CHECK(bam_extractor.ChromName(10) == "chr11");
  CHECK(bam_extractor.ChromName(2841) == "HLA-A*01:01:01:01");
}

TEST_CASE("Extractor::SetRegionBatchToExtract(absl::Span<const Reference::Region>)",
          "[lancet][hts][Extractor]") {
  Reference const ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));
  auto const case_cram_path = MakePath(FULL_DATA_DIR, CASE_CRAM_NAME);
  auto const case_bam_path = MakePath(FULL_DATA_DIR, CASE_BAM_NAME);

  auto const regions = {ref.MakeRegion("chr4:100000000-100000000"),
                        ref.MakeRegion("chr4:100000-100000")};

  SECTION("Can extract cram alignments from chr4:100000000-100000000, chr4:100000-100000") {
    Extractor cram_extractor(case_cram_path, ref, Alignment::Fields::CORE_QNAME);
    cram_extractor.SetRegionBatchToExtract(regions);
    CHECK(std::distance(cram_extractor.begin(), cram_extractor.end()) == 170);
  }

  SECTION("Can extract bam alignments from chr4:100000000-100000000, chr4:100000-100000") {
    Extractor bam_extractor(case_bam_path, ref, Alignment::Fields::CORE_QNAME);
    bam_extractor.SetRegionBatchToExtract(regions);
    CHECK(std::distance(bam_extractor.begin(), bam_extractor.end()) == 170);
  }
}

// Catch2 SECTION fan-out inflates clang-tidy's cognitive-complexity metric beyond the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Extractor::SetFilterExpression(const std::string&)", "[lancet][hts][Extractor]") {
  Reference const ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));
  auto const case_cram_path = MakePath(FULL_DATA_DIR, CASE_CRAM_NAME);
  auto const case_bam_path = MakePath(FULL_DATA_DIR, CASE_BAM_NAME);

  auto const regions = {ref.MakeRegion("chr4:100000000-100000000"),
                        ref.MakeRegion("chr4:100000-100000")};

  SECTION("Can extract cram alignments from chr4:100000000-100000000, chr4:100000-100000") {
    Extractor cram_extractor(case_cram_path, ref, Alignment::Fields::CIGAR_SEQ_QUAL);
    cram_extractor.SetRegionBatchToExtract(regions);
    cram_extractor.SetFilterExpression("min(qual) >= 20");
    CHECK(std::distance(cram_extractor.begin(), cram_extractor.end()) == 28);
  }

  SECTION("Can extract bam alignments from chr4:100000000-100000000, chr4:100000-100000") {
    Extractor bam_extractor(case_bam_path, ref, Alignment::Fields::CIGAR_SEQ_QUAL);
    bam_extractor.SetRegionBatchToExtract(regions);
    bam_extractor.SetFilterExpression("min(qual) >= 20");
    CHECK(std::distance(bam_extractor.begin(), bam_extractor.end()) == 28);
  }

  SECTION("Throws when invalid filter expression is provided to cram") {
    Extractor cram_extractor(case_cram_path, ref);
    cram_extractor.SetRegionBatchToExtract(regions);
    cram_extractor.SetFilterExpression("invalid");
    // CHECK_THROWS_AS asserts the call throws; std::distance's return value is irrelevant here.
    // NOLINTNEXTLINE(bugprone-unused-return-value)
    CHECK_THROWS_AS(std::distance(cram_extractor.begin(), cram_extractor.end()),
                    std::runtime_error);
  }

  SECTION("Throws when invalid filter expression is provided to bam") {
    Extractor bam_extractor(case_bam_path, ref);
    bam_extractor.SetRegionBatchToExtract(regions);
    bam_extractor.SetFilterExpression("invalid");
    // CHECK_THROWS_AS asserts the call throws; std::distance's return value is irrelevant here.
    // NOLINTNEXTLINE(bugprone-unused-return-value)
    CHECK_THROWS_AS(std::distance(bam_extractor.begin(), bam_extractor.end()), std::runtime_error);
  }
}
