#include "lancet/hts/alignment.h"

#include <algorithm>
#include <iterator>
#include <vector>

#include "catch_amalgamated.hpp"
#include "lancet/base/types.h"
#include "lancet/hts/extractor.h"
#include "lancet/hts/reference.h"
#include "lancet_test_config.h"

using namespace lancet::hts;

TEST_CASE("Can populate only the requested fields in bam/cram", "[lancet][hts][Alignment]") {
  const Reference ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));
  const auto tumor_cram_path = MakePath(FULL_DATA_DIR, TUMOR_CRAM_NAME);
  const auto tumor_bam_path = MakePath(FULL_DATA_DIR, TUMOR_BAM_NAME);
  const auto regions = {ref.MakeRegion("chr4:100000000")};

  SECTION("CORE_QNAME alignment fields only from cram") {
    Extractor cram_extractor(tumor_cram_path, ref, Alignment::Fields::CORE_QNAME);
    cram_extractor.SetRegionBatchToExtract(regions);
    const auto aln = *(std::begin(cram_extractor));
    CHECK_FALSE(aln.IsEmpty());
    CHECK(aln.SeqView().empty());
    CHECK(aln.QualView().empty());
    CHECK(aln.CigarData().empty());
    CHECK(aln.NumTags() == 0);
  }

  SECTION("CORE_QNAME alignment fields only from bam") {
    Extractor bam_extractor(tumor_bam_path, ref, Alignment::Fields::CORE_QNAME);
    bam_extractor.SetRegionBatchToExtract(regions);
    const auto aln = *(std::begin(bam_extractor));
    CHECK_FALSE(aln.IsEmpty());
    CHECK(aln.SeqView().empty());
    CHECK(aln.QualView().empty());
    CHECK(aln.CigarData().empty());
    CHECK(aln.NumTags() == 0);
  }

  SECTION("CIGAR_SEQ_QUAL alignment fields only from cram") {
    Extractor cram_extractor(tumor_cram_path, ref, Alignment::Fields::CIGAR_SEQ_QUAL);
    cram_extractor.SetRegionBatchToExtract(regions);
    const auto aln = *(std::begin(cram_extractor));
    CHECK_FALSE(aln.IsEmpty());
    CHECK_FALSE(aln.SeqView().empty());
    CHECK_FALSE(aln.QualView().empty());
    CHECK_FALSE(aln.CigarData().empty());
    CHECK(aln.NumTags() == 0);
  }

  SECTION("CIGAR_SEQ_QUAL alignment fields only from bam") {
    Extractor bam_extractor(tumor_bam_path, ref, Alignment::Fields::CIGAR_SEQ_QUAL);
    bam_extractor.SetRegionBatchToExtract(regions);
    const auto aln = *(std::begin(bam_extractor));
    CHECK_FALSE(aln.IsEmpty());
    CHECK_FALSE(aln.SeqView().empty());
    CHECK_FALSE(aln.QualView().empty());
    CHECK_FALSE(aln.CigarData().empty());
    CHECK(aln.NumTags() == 0);
  }

  SECTION("AUX_RGAUX alignment fields only from cram") {
    Extractor cram_extractor(tumor_cram_path, ref, Alignment::Fields::AUX_RGAUX);
    cram_extractor.SetRegionBatchToExtract(regions);
    const auto aln = *(std::begin(cram_extractor));
    CHECK_FALSE(aln.IsEmpty());
    CHECK_FALSE(aln.SeqView().empty());
    CHECK_FALSE(aln.QualView().empty());
    CHECK_FALSE(aln.CigarData().empty());
    CHECK(aln.NumTags() == 0);
  }

  SECTION("AUX_RGAUX alignment fields only from bam") {
    Extractor bam_extractor(tumor_bam_path, ref, Alignment::Fields::AUX_RGAUX);
    bam_extractor.SetRegionBatchToExtract(regions);
    const auto aln = *(std::begin(bam_extractor));
    CHECK_FALSE(aln.IsEmpty());
    CHECK_FALSE(aln.SeqView().empty());
    CHECK_FALSE(aln.QualView().empty());
    CHECK_FALSE(aln.CigarData().empty());
    CHECK(aln.NumTags() == 0);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Alignment has expected data when reading bam/cram", "[lancet][hts][Alignment]") {
  const Reference ref(MakePath(FULL_DATA_DIR, GRCH38_REF_NAME));
  const auto tumor_cram_path = MakePath(FULL_DATA_DIR, TUMOR_CRAM_NAME);
  const auto tumor_bam_path = MakePath(FULL_DATA_DIR, TUMOR_BAM_NAME);

  static constexpr auto expected_seq =
      "GAATGGAAAGGAATGGAATGGAATGGAATGGAATGGAATGGAATCAACTCGATTGCAATCGAATGGAATGGAATGGAATTAACCCGAATAGAATGGAATGGAATGGAATGGA"
      "ACGGAACGGAACGCAATGGAATGCATTGGAATGGAATGG";

  static const std::vector<u8> expected_qual{
      31, 34, 35, 34, 34, 32, 34, 34, 34, 33, 32, 34, 34, 33, 32, 32, 34, 34, 33, 32, 32, 34, 15, 33, 32, 32,
      34, 13, 33, 32, 32, 34, 34, 33, 32, 32, 34, 33, 33, 32, 32, 34, 12, 33, 33, 33, 34, 31, 34, 33, 25, 34,
      33, 35, 32, 33, 33, 34, 33, 33, 25, 34, 34, 33, 32, 32, 34, 34, 33, 32, 32, 34, 34, 33, 32, 32, 34, 34,
      33, 35, 32, 34, 31, 33, 33, 25, 34, 34, 33, 32, 33, 34, 13, 33, 32, 32, 34, 13, 33, 32, 32, 35, 35, 33,
      32, 32, 34, 35, 33, 32, 32, 35, 13, 31, 25, 32, 34, 13, 31, 26, 33, 35, 13, 31, 26, 34, 34, 13, 34, 33,
      32, 35, 35, 34, 33, 34, 11, 34, 36, 33, 33, 35, 26, 34, 34, 34, 26, 26, 35, 33, 11};

  static const std::vector<std::string_view> expected_tags{"AS", "MC", "MD", "NM", "RG", "SA", "XS", "pa"};

  SECTION("Can read cram alignment records") {
    Extractor cram_extractor(tumor_cram_path, ref, Alignment::Fields::AUX_RGAUX,
                             {"RG", "MC", "NM", "SA", "XS", "MD", "AS", "pa"});

    cram_extractor.SetFilterExpression("mapq >= 20 && tlen >= 300");
    const auto aln = *(std::begin(cram_extractor));
    CHECK(aln.StartPos0() == 789838);
    CHECK(aln.MateStartPos0() == 41877838);
    CHECK(aln.InsertSize() == 7318);
    CHECK(aln.ChromIndex() == 0);
    CHECK(aln.MateChromIndex() == 9);
    CHECK(aln.Flag().HasFlagsSet(2145));
    CHECK(aln.MapQual() == 21);
    CHECK(aln.QnameView() == "SRR7890893.726153849");
    CHECK(aln.SeqView() == expected_seq);
    CHECK(aln.QualView() == expected_qual);
    CHECK(aln.CigarString() == "62S33M5D39M17S");
    CHECK(aln.TagNamesView() == expected_tags);
    std::for_each(expected_tags.cbegin(), expected_tags.cend(),
                  [&aln](const auto& tag_name) { CHECK(aln.HasTag(tag_name)); });

    CHECK(aln.GetTag<i64>("AS").ok());
    CHECK(aln.GetTag<i64>("AS").value() == 31);

    CHECK(aln.GetTag<std::string_view>("MC").ok());
    CHECK(aln.GetTag<std::string_view>("MC").value() == "22S38M91S");

    CHECK(aln.GetTag<std::string_view>("MD").ok());
    CHECK(aln.GetTag<std::string_view>("MD").value() == "17G0G2T5G5^GAACC3A9C25");

    CHECK(aln.GetTag<i64>("NM").ok());
    CHECK(aln.GetTag<i64>("NM").value() == 11);

    CHECK(aln.GetTag<std::string_view>("RG").ok());
    CHECK(aln.GetTag<std::string_view>("RG").value() == "SRR7890893rg");

    CHECK(aln.GetTag<std::string_view>("SA").ok());
    CHECK(aln.GetTag<std::string_view>("SA").value() ==
          "chr10,41870559,+,49M102S,60,0;chrUn_KN707896v1_decoy,13830,-,151M,19,4;");

    CHECK(aln.GetTag<i64>("XS").ok());
    CHECK(aln.GetTag<i64>("XS").value() == 0);

    CHECK(aln.GetTag<f64>("pa").ok());
    static constexpr double EXPECTED_PA = 0.237;
    static constexpr double epsilon = 0.001;
    CHECK_THAT(aln.GetTag<f64>("pa").value(), Catch::Matchers::WithinRel(EXPECTED_PA, epsilon));
  }

  SECTION("Can read bam alignment records") {
    Extractor bam_extractor(tumor_bam_path, ref, Alignment::Fields::AUX_RGAUX,
                            {"RG", "MC", "NM", "SA", "XS", "MD", "AS", "pa"});
    bam_extractor.SetFilterExpression("mapq >= 20 && tlen >= 300");
    const auto aln = *(std::begin(bam_extractor));
    CHECK(aln.StartPos0() == 789838);
    CHECK(aln.MateStartPos0() == 41877838);
    CHECK(aln.InsertSize() == 7318);
    CHECK(aln.ChromIndex() == 0);
    CHECK(aln.MateChromIndex() == 9);
    CHECK(aln.Flag().HasFlagsSet(2145));
    CHECK(aln.MapQual() == 21);
    CHECK(aln.QnameView() == "SRR7890893.726153849");
    CHECK(aln.SeqView() == expected_seq);
    CHECK(aln.QualView() == expected_qual);
    CHECK(aln.CigarString() == "62S33M5D39M17S");
    CHECK(aln.TagNamesView() == expected_tags);
    std::for_each(expected_tags.cbegin(), expected_tags.cend(),
                  [&aln](const auto& tag_name) { CHECK(aln.HasTag(tag_name)); });

    CHECK(aln.GetTag<i64>("AS").ok());
    CHECK(aln.GetTag<i64>("AS").value() == 31);

    CHECK(aln.GetTag<std::string_view>("MC").ok());
    CHECK(aln.GetTag<std::string_view>("MC").value() == "22S38M91S");

    CHECK(aln.GetTag<std::string_view>("MD").ok());
    CHECK(aln.GetTag<std::string_view>("MD").value() == "17G0G2T5G5^GAACC3A9C25");

    CHECK(aln.GetTag<i64>("NM").ok());
    CHECK(aln.GetTag<i64>("NM").value() == 11);

    CHECK(aln.GetTag<std::string_view>("RG").ok());
    CHECK(aln.GetTag<std::string_view>("RG").value() == "SRR7890893rg");

    CHECK(aln.GetTag<std::string_view>("SA").ok());
    CHECK(aln.GetTag<std::string_view>("SA").value() ==
          "chr10,41870559,+,49M102S,60,0;chrUn_KN707896v1_decoy,13830,-,151M,19,4;");

    CHECK(aln.GetTag<i64>("XS").ok());
    CHECK(aln.GetTag<i64>("XS").value() == 0);

    CHECK(aln.GetTag<f64>("pa").ok());
    static constexpr double EXPECTED_PA = 0.237;
    static constexpr double epsilon = 0.001;
    CHECK_THAT(aln.GetTag<f64>("pa").value(), Catch::Matchers::WithinRel(EXPECTED_PA, epsilon));
  }
}
