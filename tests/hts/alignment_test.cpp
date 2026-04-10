#include "lancet/hts/alignment.h"

#include "lancet/base/types.h"
#include "lancet/hts/extractor.h"
#include "lancet/hts/reference.h"

#include "catch_amalgamated.hpp"
#include "lancet_test_config.h"

#include <algorithm>
#include <iterator>
#include <vector>

using lancet::hts::Alignment;
using lancet::hts::Extractor;
using lancet::hts::Reference;

TEST_CASE("Can populate only the requested fields in bam/cram", "[lancet][hts][Alignment]") {
  Reference const ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));
  auto const tumor_cram_path = MakePath(FULL_DATA_DIR, TUMOR_CRAM_NAME);
  auto const tumor_bam_path = MakePath(FULL_DATA_DIR, TUMOR_BAM_NAME);
  auto const regions = {ref.MakeRegion("chr4:100000000")};

  // NOTE: With the zero-copy Alignment proxy, BuildSequence/BuildQualities/CigarData
  // may return non-empty data even when only CORE_QNAME is requested, because htslib
  // may populate the full bam1_t record in memory (especially for BAM). The Fields
  // enum controls what htslib *guarantees* to populate, not what it excludes.
  SECTION("CORE_QNAME alignment fields only from cram") {
    Extractor cram_extractor(tumor_cram_path, ref, Alignment::Fields::CORE_QNAME);
    cram_extractor.SetRegionBatchToExtract(regions);
    auto const align = *(std::begin(cram_extractor));
    CHECK_FALSE(align.IsEmpty());
    CHECK_FALSE(align.QnameView().empty());
  }

  SECTION("CORE_QNAME alignment fields only from bam") {
    Extractor bam_extractor(tumor_bam_path, ref, Alignment::Fields::CORE_QNAME);
    bam_extractor.SetRegionBatchToExtract(regions);
    auto const align = *(std::begin(bam_extractor));
    CHECK_FALSE(align.IsEmpty());
    CHECK_FALSE(align.QnameView().empty());
  }

  SECTION("CIGAR_SEQ_QUAL alignment fields only from cram") {
    Extractor cram_extractor(tumor_cram_path, ref, Alignment::Fields::CIGAR_SEQ_QUAL);
    cram_extractor.SetRegionBatchToExtract(regions);
    auto const align = *(std::begin(cram_extractor));
    CHECK_FALSE(align.IsEmpty());
    CHECK_FALSE(align.BuildSequence().empty());
    CHECK_FALSE(align.BuildQualities().empty());
    CHECK_FALSE(align.CigarData().empty());
  }

  SECTION("CIGAR_SEQ_QUAL alignment fields only from bam") {
    Extractor bam_extractor(tumor_bam_path, ref, Alignment::Fields::CIGAR_SEQ_QUAL);
    bam_extractor.SetRegionBatchToExtract(regions);
    auto const align = *(std::begin(bam_extractor));
    CHECK_FALSE(align.IsEmpty());
    CHECK_FALSE(align.BuildSequence().empty());
    CHECK_FALSE(align.BuildQualities().empty());
    CHECK_FALSE(align.CigarData().empty());
  }

  SECTION("AUX_RGAUX alignment fields only from cram") {
    Extractor cram_extractor(tumor_cram_path, ref, Alignment::Fields::AUX_RGAUX);
    cram_extractor.SetRegionBatchToExtract(regions);
    auto const align = *(std::begin(cram_extractor));
    CHECK_FALSE(align.IsEmpty());
    CHECK_FALSE(align.BuildSequence().empty());
    CHECK_FALSE(align.BuildQualities().empty());
    CHECK_FALSE(align.CigarData().empty());
  }

  SECTION("AUX_RGAUX alignment fields only from bam") {
    Extractor bam_extractor(tumor_bam_path, ref, Alignment::Fields::AUX_RGAUX);
    bam_extractor.SetRegionBatchToExtract(regions);
    auto const align = *(std::begin(bam_extractor));
    CHECK_FALSE(align.IsEmpty());
    CHECK_FALSE(align.BuildSequence().empty());
    CHECK_FALSE(align.BuildQualities().empty());
    CHECK_FALSE(align.CigarData().empty());
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Alignment has expected data when reading bam/cram", "[lancet][hts][Alignment]") {
  Reference const ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));
  auto const tumor_cram_path = MakePath(FULL_DATA_DIR, TUMOR_CRAM_NAME);
  auto const tumor_bam_path = MakePath(FULL_DATA_DIR, TUMOR_BAM_NAME);

  static constexpr auto EXPECTED_SEQ = "GAATGGAAAGGAATGGAATGGAATGGAATGGAATGGAATGGAATCAACTCGATTGCAAT"
                                       "CGAATGGAATGGAATGGAATTAACCCGAATAGAATGGAATGGAATGGAATGGA"
                                       "ACGGAACGGAACGCAATGGAATGCATTGGAATGGAATGG";

  static std::vector<u8> const EXPECTED_QUAL{
      31, 34, 35, 34, 34, 32, 34, 34, 34, 33, 32, 34, 34, 33, 32, 32, 34, 34, 33, 32, 32, 34,
      15, 33, 32, 32, 34, 13, 33, 32, 32, 34, 34, 33, 32, 32, 34, 33, 33, 32, 32, 34, 12, 33,
      33, 33, 34, 31, 34, 33, 25, 34, 33, 35, 32, 33, 33, 34, 33, 33, 25, 34, 34, 33, 32, 32,
      34, 34, 33, 32, 32, 34, 34, 33, 32, 32, 34, 34, 33, 35, 32, 34, 31, 33, 33, 25, 34, 34,
      33, 32, 33, 34, 13, 33, 32, 32, 34, 13, 33, 32, 32, 35, 35, 33, 32, 32, 34, 35, 33, 32,
      32, 35, 13, 31, 25, 32, 34, 13, 31, 26, 33, 35, 13, 31, 26, 34, 34, 13, 34, 33, 32, 35,
      35, 34, 33, 34, 11, 34, 36, 33, 33, 35, 26, 34, 34, 34, 26, 26, 35, 33, 11};

  static std::vector<std::string_view> const EXPECTED_TAGS{"AS", "MC", "MD", "NM",
                                                           "RG", "SA", "XS", "pa"};

  SECTION("Can read cram alignment records") {
    Extractor cram_extractor(tumor_cram_path, ref, Alignment::Fields::AUX_RGAUX,
                             {"RG", "MC", "NM", "SA", "XS", "MD", "AS", "pa"});

    cram_extractor.SetFilterExpression("mapq >= 20 && tlen >= 300");
    auto const align = *(std::begin(cram_extractor));
    CHECK(align.StartPos0() == 789'838);
    CHECK(align.MateStartPos0() == 41'877'838);
    CHECK(align.InsertSize() == 7318);
    CHECK(align.ChromIndex() == 0);
    CHECK(align.MateChromIndex() == 9);
    CHECK(align.Flag().HasFlagsSet(2145));
    CHECK(align.MapQual() == 21);
    CHECK(align.QnameView() == "SRR7890893.726153849");
    CHECK(align.BuildSequence() == EXPECTED_SEQ);
    CHECK(align.BuildQualities() == EXPECTED_QUAL);
    CHECK(align.CigarString() == "62S33M5D39M17S");
    std::for_each(EXPECTED_TAGS.cbegin(), EXPECTED_TAGS.cend(),
                  [&align](auto const& tag_name) { CHECK(align.HasTag(tag_name)); });

    CHECK(align.GetTag<i64>("AS").ok());
    CHECK(align.GetTag<i64>("AS").value() == 31);

    CHECK(align.GetTag<std::string_view>("MC").ok());
    CHECK(align.GetTag<std::string_view>("MC").value() == "22S38M91S");

    CHECK(align.GetTag<std::string_view>("MD").ok());
    CHECK(align.GetTag<std::string_view>("MD").value() == "17G0G2T5G5^GAACC3A9C25");

    CHECK(align.GetTag<i64>("NM").ok());
    CHECK(align.GetTag<i64>("NM").value() == 11);

    CHECK(align.GetTag<std::string_view>("RG").ok());
    CHECK(align.GetTag<std::string_view>("RG").value() == "SRR7890893rg");

    CHECK(align.GetTag<std::string_view>("SA").ok());
    CHECK(align.GetTag<std::string_view>("SA").value() ==
          "chr10,41870559,+,49M102S,60,0;chrUn_KN707896v1_decoy,13830,-,151M,19,4;");

    CHECK(align.GetTag<i64>("XS").ok());
    CHECK(align.GetTag<i64>("XS").value() == 0);

    CHECK(align.GetTag<f64>("pa").ok());
    static constexpr double EXPECTED_PA = 0.237;
    static constexpr double EPSILON = 0.001;
    CHECK_THAT(align.GetTag<f64>("pa").value(), Catch::Matchers::WithinRel(EXPECTED_PA, EPSILON));
  }

  SECTION("Can read bam alignment records") {
    Extractor bam_extractor(tumor_bam_path, ref, Alignment::Fields::AUX_RGAUX,
                            {"RG", "MC", "NM", "SA", "XS", "MD", "AS", "pa"});
    bam_extractor.SetFilterExpression("mapq >= 20 && tlen >= 300");
    auto const align = *(std::begin(bam_extractor));
    CHECK(align.StartPos0() == 789'838);
    CHECK(align.MateStartPos0() == 41'877'838);
    CHECK(align.InsertSize() == 7318);
    CHECK(align.ChromIndex() == 0);
    CHECK(align.MateChromIndex() == 9);
    CHECK(align.Flag().HasFlagsSet(2145));
    CHECK(align.MapQual() == 21);
    CHECK(align.QnameView() == "SRR7890893.726153849");
    CHECK(align.BuildSequence() == EXPECTED_SEQ);
    CHECK(align.BuildQualities() == EXPECTED_QUAL);
    CHECK(align.CigarString() == "62S33M5D39M17S");
    std::for_each(EXPECTED_TAGS.cbegin(), EXPECTED_TAGS.cend(),
                  [&align](auto const& tag_name) { CHECK(align.HasTag(tag_name)); });

    CHECK(align.GetTag<i64>("AS").ok());
    CHECK(align.GetTag<i64>("AS").value() == 31);

    CHECK(align.GetTag<std::string_view>("MC").ok());
    CHECK(align.GetTag<std::string_view>("MC").value() == "22S38M91S");

    CHECK(align.GetTag<std::string_view>("MD").ok());
    CHECK(align.GetTag<std::string_view>("MD").value() == "17G0G2T5G5^GAACC3A9C25");

    CHECK(align.GetTag<i64>("NM").ok());
    CHECK(align.GetTag<i64>("NM").value() == 11);

    CHECK(align.GetTag<std::string_view>("RG").ok());
    CHECK(align.GetTag<std::string_view>("RG").value() == "SRR7890893rg");

    CHECK(align.GetTag<std::string_view>("SA").ok());
    CHECK(align.GetTag<std::string_view>("SA").value() ==
          "chr10,41870559,+,49M102S,60,0;chrUn_KN707896v1_decoy,13830,-,151M,19,4;");

    CHECK(align.GetTag<i64>("XS").ok());
    CHECK(align.GetTag<i64>("XS").value() == 0);

    CHECK(align.GetTag<f64>("pa").ok());
    static constexpr double EXPECTED_PA = 0.237;
    static constexpr double EPSILON = 0.001;
    CHECK_THAT(align.GetTag<f64>("pa").value(), Catch::Matchers::WithinRel(EXPECTED_PA, EPSILON));
  }
}
