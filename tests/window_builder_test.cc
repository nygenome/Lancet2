#include "lancet2/window_builder.h"

#include "absl/strings/str_format.h"
#include "catch2/catch_test_macros.hpp"
#include "generated/test_config.h"

TEST_CASE("Can build windows flexibly from various sources", "[lancet2::WindowBuilder]") {
  const auto refPath = absl::StrFormat("%s/human_g1k_v37.1_1_90000000.fa.gz", TEST_DATA_DIR);
  const auto bedPath = absl::StrFormat("%s/example.bed", TEST_DATA_DIR);
  lancet2::WindowBuilder wb(refPath, 250, 500, 50);

  CHECK(wb.GetSize() == 0);
  CHECK(wb.IsEmpty());
  CHECK(lancet2::WindowBuilder::StepSize(84, 500) == 100);

  CHECK(wb.AddRegion("1:1-3000").ok());
  CHECK(wb.AddRegion("1:1-3000").ok());
  CHECK(wb.AddRegion("1:1-3000").ok());
  CHECK(wb.AddRegionsBatch({"1:1-3000", "1:1-3000"}).ok());
  CHECK(wb.AddRegionsFromBedFile(bedPath).ok());
  const auto r1 = wb.BuildWindows();
  CHECK(r1.ok());
  CHECK(r1.value().size() == 11);

  wb.AddAllRefRegions();
  const auto r2 = wb.BuildWindows();
  CHECK(r2.ok());
  CHECK(r2.value().size() == 300000);
}
