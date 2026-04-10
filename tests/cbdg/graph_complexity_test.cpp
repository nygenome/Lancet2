#include "lancet/cbdg/graph_complexity.h"

#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"

#include "absl/types/span.h"
#include "catch_amalgamated.hpp"

#include <vector>

#include <cmath>

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  GraphComplexity::IsComplex Thresholds                                   ║
// ╚══════════════════════════════════════════════════════════════════════════╝

// NOTE: We cannot directly construct GraphComplexity with specific field values
// because the fields are private (populated only by ComputeGraphComplexity).
// IsComplex() tests use the default-constructed state plus the public constants.

TEST_CASE("GraphComplexity::IsComplex default is not complex", "[lancet][cbdg][graph_complexity]") {
  using lancet::cbdg::GraphComplexity;

  GraphComplexity const cplx;
  CHECK_FALSE(cplx.IsComplex());
  CHECK(cplx.CyclomaticComplexity() == 0);
  CHECK(cplx.NumBranchPoints() == 0);
}

TEST_CASE("GraphComplexity::IsComplex thresholds are correct", "[lancet][cbdg][graph_complexity]") {
  using lancet::cbdg::GraphComplexity;

  CHECK(GraphComplexity::MAX_CYCLOMATIC_COMPLEXITY == 50);
  CHECK(GraphComplexity::MAX_BRANCH_POINTS == 50);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  GraphComplexity::GraphEntanglementIndex (GEI)                           ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("GEI: default-constructed → GEI = 0", "[lancet][cbdg][graph_complexity][gei]") {
  using lancet::cbdg::GraphComplexity;

  GraphComplexity const cplx;
  CHECK(cplx.GraphEntanglementIndex() == Catch::Approx(0.0).margin(1e-10));
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  GraphMetrics::FormatVcfValue                                            ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("GraphMetrics::FormatVcfValue: 3 comma-separated values",
          "[lancet][cbdg][graph_complexity][format]") {
  lancet::caller::RawVariant::GraphMetrics gmtx;
  gmtx.mGraphEntanglementIndex = 2.345;
  gmtx.mTipToPathCovRatio = 0.15;
  gmtx.mMaxSingleDirDegree = 4;

  auto const vcf_val = gmtx.FormatVcfValue();
  // 3 values → 2 commas
  CHECK(std::count(vcf_val.begin(), vcf_val.end(), ',') == 2);
  // Ends with ",4" (the MaxSingleDirDegree)
  CHECK(vcf_val.find(",4") != std::string::npos);
}
