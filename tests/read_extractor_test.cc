#include "lancet2/read_extractor.h"

#include "absl/strings/str_format.h"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "generated/test_config.h"

TEST_CASE("Can extract reads from region in a bam file", "[lancet2::ReadExtractor]") {
  const auto bamPath = absl::StrFormat("%s/human_g1k_v37.no_mutations.1_82960000_82970000.bam", TEST_DATA_DIR);
  const auto refPath = absl::StrFormat("%s/human_g1k_v37.1_1_90000000.fa.gz", TEST_DATA_DIR);
  auto params = std::make_shared<lancet2::CliParams>();
  params->tumorPath = bamPath;
  params->normalPath = bamPath;
  params->referencePath = refPath;

  lancet2::ReadExtractor re(std::const_pointer_cast<const lancet2::CliParams>(params));
  const auto scanResult = re.ScanRegion(lancet2::GenomicRegion{"1", 82960000, 82970000});
  CHECK_FALSE(scanResult.HasMutationEvidence);
  CHECK_THAT(scanResult.AverageCoverage, Catch::Matchers::WithinRel(200, 0.01));
  const auto reads = re.ExtractReads(lancet2::GenomicRegion{"1", 82960000, 82970000});
  CHECK(reads.size() == 12830);
}

TEST_CASE("Can extract reads from region in a cram file", "[lancet2::ReadExtractor]") {
  const auto tumorCram = absl::StrFormat("%s/tumor.cram", TEST_DATA_DIR);
  const auto normalCram = absl::StrFormat("%s/normal.cram", TEST_DATA_DIR);
  const auto refPath = absl::StrFormat("%s/human_g1k_v37.1_1_90000000.fa.gz", TEST_DATA_DIR);

  auto params = std::make_shared<lancet2::CliParams>();
  params->tumorPath = tumorCram;
  params->normalPath = normalCram;
  params->referencePath = refPath;

  lancet2::ReadExtractor re(std::const_pointer_cast<const lancet2::CliParams>(params));
  const auto scanResult = re.ScanRegion(lancet2::GenomicRegion{"1", 82960000, 82970000});
  CHECK(scanResult.HasMutationEvidence);
  CHECK_THAT(scanResult.AverageCoverage, Catch::Matchers::WithinRel(151, 0.01));
  const auto reads = re.ExtractReads(lancet2::GenomicRegion{"1", 82960000, 82970000});
  CHECK(reads.size() == 9529);
}
