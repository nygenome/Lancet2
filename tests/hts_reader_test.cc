#include "lancet2/hts_reader.h"

#include "absl/strings/str_format.h"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_vector.hpp"
#include "generated/test_config.h"

struct AlnFileStats {
  int Duplicate = 0;
  int Supplementary = 0;
  int Primary = 0;
  int Secondary = 0;
  int QcFailed = 0;
  int Unmapped = 0;
  int MateUnmapped = 0;
  int ReverseStrand = 0;
  int ForwardStrand = 0;
  int MateReverseStrand = 0;
  int MateForwardStrand = 0;
  int Paired = 0;
  int ProperPaired = 0;
  int Read1 = 0;
  int Read2 = 0;
};

TEST_CASE("Can read BAM files", "[lancet::hts::HtsReader]") {
  const auto bamPath = absl::StrFormat("%s/tumor.bam", TEST_DATA_DIR);
  const auto refPath = absl::StrFormat("%s/human_g1k_v37.1_1_90000000.fa.gz", TEST_DATA_DIR);
  lancet2::HtsReader rdr(bamPath, refPath);

  SECTION("can get sample name from bam header") { CHECK(rdr.GetSampleName() == "NA12892_mem_binomial_indel"); }
  SECTION("can get contig index from bam header") {
    CHECK(rdr.GetContigIndex("1") == 0);
    CHECK(rdr.GetContigIndex("10") == 9);
  }

  SECTION("can check for existence of bam tags") {
    CHECK(lancet2::TagPeekCheck(bamPath, refPath, "MD"));
    CHECK_FALSE(lancet2::TagPeekCheck(bamPath, refPath, "BX"));
  }

  SECTION("can set contig and jump with iterator") { CHECK(rdr.JumpToContig("1").ok()); }
  SECTION("can set region and jump with iterator") {
    CHECK(rdr.JumpToRegion(lancet2::GenomicRegion{"1", 82960000, 82970000}).ok());
  }
  SECTION("can reset iterator to begining of file") { CHECK_NOTHROW(rdr.ResetIterator()); }

  SECTION("aggregate read statistics from bam match expected values") {
    rdr.ResetIterator();
    lancet2::HtsAlignment aln;
    AlnFileStats stats;
    auto state = rdr.GetNextAlignment(&aln, {});
    while (state == lancet2::HtsReader::IteratorState::VALID) {
      if (aln.IsDuplicate()) stats.Duplicate += 1;
      if (aln.IsSupplementary()) stats.Supplementary += 1;
      if (aln.IsPrimary()) stats.Primary += 1;
      if (aln.IsSecondary()) stats.Secondary += 1;
      if (aln.IsQcFailed()) stats.QcFailed += 1;
      if (aln.IsUnmapped()) stats.Unmapped += 1;
      if (aln.IsMateUnmapped()) stats.MateUnmapped += 1;
      if (aln.IsReverseStrand()) stats.ReverseStrand += 1;
      if (!aln.IsReverseStrand()) stats.ForwardStrand += 1;
      if (aln.IsMateReverseStrand()) stats.MateReverseStrand += 1;
      if (!aln.IsMateReverseStrand()) stats.MateForwardStrand += 1;
      if (aln.IsPaired()) stats.Paired += 1;
      if (aln.IsProperPair()) stats.ProperPaired += 1;
      if (aln.IsRead1()) stats.Read1 += 1;
      if (aln.IsRead2()) stats.Read2 += 2;

      state = rdr.GetNextAlignment(&aln, {});
    }

    CHECK(state == lancet2::HtsReader::IteratorState::DONE);
    CHECK(stats.Duplicate == 563);
    CHECK(stats.Supplementary == 0);
    CHECK(stats.Primary == 7187);
    CHECK(stats.Secondary == 0);
    CHECK(stats.QcFailed == 0);
    CHECK(stats.Unmapped == 11);
    CHECK(stats.MateUnmapped == 11);
    CHECK(stats.ReverseStrand == 3565);
    CHECK(stats.ForwardStrand == 3622);
    CHECK(stats.MateReverseStrand == 3593);
    CHECK(stats.MateForwardStrand == 3594);
    CHECK(stats.Paired == 7187);
    CHECK(stats.ProperPaired == 7100);
    CHECK(stats.Read1 == 3602);
    CHECK(stats.Read2 == 7170);
  }

  SECTION("can get contig infos from bam header") {
    CHECK_THAT(rdr.GetContigs(),
               Catch::Matchers::Equals(std::vector<lancet2::ContigInfo>{
                   {"1", 249250621},       {"2", 243199373},       {"3", 198022430},       {"4", 191154276},
                   {"5", 180915260},       {"6", 171115067},       {"7", 159138663},       {"8", 146364022},
                   {"9", 141213431},       {"10", 135534747},      {"11", 135006516},      {"12", 133851895},
                   {"13", 115169878},      {"14", 107349540},      {"15", 102531392},      {"16", 90354753},
                   {"17", 81195210},       {"18", 78077248},       {"19", 59128983},       {"20", 63025520},
                   {"21", 48129895},       {"22", 51304566},       {"X", 155270560},       {"Y", 59373566},
                   {"MT", 16569},          {"GL000207.1", 4262},   {"GL000226.1", 15008},  {"GL000229.1", 19913},
                   {"GL000231.1", 27386},  {"GL000210.1", 27682},  {"GL000239.1", 33824},  {"GL000235.1", 34474},
                   {"GL000201.1", 36148},  {"GL000247.1", 36422},  {"GL000245.1", 36651},  {"GL000197.1", 37175},
                   {"GL000203.1", 37498},  {"GL000246.1", 38154},  {"GL000249.1", 38502},  {"GL000196.1", 38914},
                   {"GL000248.1", 39786},  {"GL000244.1", 39929},  {"GL000238.1", 39939},  {"GL000202.1", 40103},
                   {"GL000234.1", 40531},  {"GL000232.1", 40652},  {"GL000206.1", 41001},  {"GL000240.1", 41933},
                   {"GL000236.1", 41934},  {"GL000241.1", 42152},  {"GL000243.1", 43341},  {"GL000242.1", 43523},
                   {"GL000230.1", 43691},  {"GL000237.1", 45867},  {"GL000233.1", 45941},  {"GL000204.1", 81310},
                   {"GL000198.1", 90085},  {"GL000208.1", 92689},  {"GL000191.1", 106433}, {"GL000227.1", 128374},
                   {"GL000228.1", 129120}, {"GL000214.1", 137718}, {"GL000221.1", 155397}, {"GL000209.1", 159169},
                   {"GL000218.1", 161147}, {"GL000220.1", 161802}, {"GL000213.1", 164239}, {"GL000211.1", 166566},
                   {"GL000199.1", 169874}, {"GL000217.1", 172149}, {"GL000216.1", 172294}, {"GL000215.1", 172545},
                   {"GL000205.1", 174588}, {"GL000219.1", 179198}, {"GL000224.1", 179693}, {"GL000223.1", 180455},
                   {"GL000195.1", 182896}, {"GL000212.1", 186858}, {"GL000222.1", 186861}, {"GL000200.1", 187035},
                   {"GL000193.1", 189789}, {"GL000194.1", 191469}, {"GL000225.1", 211173}, {"GL000192.1", 547496},
                   {"NC_007605", 171823},  {"hs37d5", 35477943}}));
  }
}

TEST_CASE("Can read CRAM files", "[lancet::hts::HtsReader]") {
  const auto cramPath = absl::StrFormat("%s/tumor.cram", TEST_DATA_DIR);
  const auto refPath = absl::StrFormat("%s/human_g1k_v37.1_1_90000000.fa.gz", TEST_DATA_DIR);
  lancet2::HtsReader rdr(cramPath, refPath);

  SECTION("can get sample name from cram header") { CHECK(rdr.GetSampleName() == "NA12892_mem_binomial_indel"); }
  SECTION("can get contig index from cram header") {
    CHECK(rdr.GetContigIndex("1") == 0);
    CHECK(rdr.GetContigIndex("10") == 9);
  }

  SECTION("can check for existence of cram tags") {
    CHECK(lancet2::TagPeekCheck(cramPath, refPath, "MD"));
    CHECK_FALSE(lancet2::TagPeekCheck(cramPath, refPath, "BX"));
  }

  SECTION("can set contig and jump with iterator") { CHECK(rdr.JumpToContig("1").ok()); }
  SECTION("can set region and jump with iterator") {
    CHECK(rdr.JumpToRegion(lancet2::GenomicRegion{"1", 82960000, 82970000}).ok());
  }
  SECTION("can reset iterator to begining of file") { CHECK_NOTHROW(rdr.ResetIterator()); }

  SECTION("aggregate read statistics from cram match expected values") {
    rdr.ResetIterator();
    lancet2::HtsAlignment aln;
    AlnFileStats stats;
    auto state = rdr.GetNextAlignment(&aln, {});
    while (state == lancet2::HtsReader::IteratorState::VALID) {
      if (aln.IsDuplicate()) stats.Duplicate += 1;
      if (aln.IsSupplementary()) stats.Supplementary += 1;
      if (aln.IsPrimary()) stats.Primary += 1;
      if (aln.IsSecondary()) stats.Secondary += 1;
      if (aln.IsQcFailed()) stats.QcFailed += 1;
      if (aln.IsUnmapped()) stats.Unmapped += 1;
      if (aln.IsMateUnmapped()) stats.MateUnmapped += 1;
      if (aln.IsReverseStrand()) stats.ReverseStrand += 1;
      if (!aln.IsReverseStrand()) stats.ForwardStrand += 1;
      if (aln.IsMateReverseStrand()) stats.MateReverseStrand += 1;
      if (!aln.IsMateReverseStrand()) stats.MateForwardStrand += 1;
      if (aln.IsPaired()) stats.Paired += 1;
      if (aln.IsProperPair()) stats.ProperPaired += 1;
      if (aln.IsRead1()) stats.Read1 += 1;
      if (aln.IsRead2()) stats.Read2 += 2;

      state = rdr.GetNextAlignment(&aln, {});
    }

    CHECK(state == lancet2::HtsReader::IteratorState::DONE);
    CHECK(stats.Duplicate == 563);
    CHECK(stats.Supplementary == 0);
    CHECK(stats.Primary == 7187);
    CHECK(stats.Secondary == 0);
    CHECK(stats.QcFailed == 0);
    CHECK(stats.Unmapped == 11);
    CHECK(stats.MateUnmapped == 11);
    CHECK(stats.ReverseStrand == 3565);
    CHECK(stats.ForwardStrand == 3622);
    CHECK(stats.MateReverseStrand == 3593);
    CHECK(stats.MateForwardStrand == 3594);
    CHECK(stats.Paired == 7187);
    CHECK(stats.ProperPaired == 7100);
    CHECK(stats.Read1 == 3602);
    CHECK(stats.Read2 == 7170);
  }

  SECTION("can get contig infos from cram header") {
    CHECK_THAT(rdr.GetContigs(),
               Catch::Matchers::Equals(std::vector<lancet2::ContigInfo>{
                   {"1", 90000000},        {"2", 243199373},       {"3", 198022430},       {"4", 191154276},
                   {"5", 180915260},       {"6", 171115067},       {"7", 159138663},       {"8", 146364022},
                   {"9", 141213431},       {"10", 135534747},      {"11", 135006516},      {"12", 133851895},
                   {"13", 115169878},      {"14", 107349540},      {"15", 102531392},      {"16", 90354753},
                   {"17", 81195210},       {"18", 78077248},       {"19", 59128983},       {"20", 63025520},
                   {"21", 48129895},       {"22", 51304566},       {"X", 155270560},       {"Y", 59373566},
                   {"MT", 16569},          {"GL000207.1", 4262},   {"GL000226.1", 15008},  {"GL000229.1", 19913},
                   {"GL000231.1", 27386},  {"GL000210.1", 27682},  {"GL000239.1", 33824},  {"GL000235.1", 34474},
                   {"GL000201.1", 36148},  {"GL000247.1", 36422},  {"GL000245.1", 36651},  {"GL000197.1", 37175},
                   {"GL000203.1", 37498},  {"GL000246.1", 38154},  {"GL000249.1", 38502},  {"GL000196.1", 38914},
                   {"GL000248.1", 39786},  {"GL000244.1", 39929},  {"GL000238.1", 39939},  {"GL000202.1", 40103},
                   {"GL000234.1", 40531},  {"GL000232.1", 40652},  {"GL000206.1", 41001},  {"GL000240.1", 41933},
                   {"GL000236.1", 41934},  {"GL000241.1", 42152},  {"GL000243.1", 43341},  {"GL000242.1", 43523},
                   {"GL000230.1", 43691},  {"GL000237.1", 45867},  {"GL000233.1", 45941},  {"GL000204.1", 81310},
                   {"GL000198.1", 90085},  {"GL000208.1", 92689},  {"GL000191.1", 106433}, {"GL000227.1", 128374},
                   {"GL000228.1", 129120}, {"GL000214.1", 137718}, {"GL000221.1", 155397}, {"GL000209.1", 159169},
                   {"GL000218.1", 161147}, {"GL000220.1", 161802}, {"GL000213.1", 164239}, {"GL000211.1", 166566},
                   {"GL000199.1", 169874}, {"GL000217.1", 172149}, {"GL000216.1", 172294}, {"GL000215.1", 172545},
                   {"GL000205.1", 174588}, {"GL000219.1", 179198}, {"GL000224.1", 179693}, {"GL000223.1", 180455},
                   {"GL000195.1", 182896}, {"GL000212.1", 186858}, {"GL000222.1", 186861}, {"GL000200.1", 187035},
                   {"GL000193.1", 189789}, {"GL000194.1", 191469}, {"GL000225.1", 211173}, {"GL000192.1", 547496},
                   {"NC_007605", 171823},  {"hs37d5", 35477943}}));
  }
}
