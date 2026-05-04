#include "lancet/cbdg/graph_complexity.h"

#include "lancet/caller/raw_variant.h"

#include "catch_amalgamated.hpp"

#include <algorithm>
#include <string>

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  GraphComplexity::IsComplex Thresholds                                   ║
// ╚══════════════════════════════════════════════════════════════════════════╝

// NOTE: We cannot directly construct GraphComplexity with specific field values
// because the fields are private (populated only by ComputeGraphComplexity).
// IsComplex() tests use the default-constructed state plus the public constants.

TEST_CASE("GraphComplexity::IsComplex default is not complex", "[lancet][cbdg][GraphComplexity]") {
  using lancet::cbdg::GraphComplexity;

  GraphComplexity const cplx;
  CHECK_FALSE(cplx.IsComplex());
  CHECK(cplx.CyclomaticComplexity() == 0);
  CHECK(cplx.NumBranchPoints() == 0);
}

TEST_CASE("GraphComplexity::IsComplex thresholds are correct", "[lancet][cbdg][GraphComplexity]") {
  using lancet::cbdg::GraphComplexity;

  CHECK(GraphComplexity::MAX_CYCLOMATIC_COMPLEXITY == 50);
  CHECK(GraphComplexity::MAX_BRANCH_POINTS == 50);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  GraphComplexity::GraphEntanglementIndex (GEI)                           ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("GEI: default-constructed → GEI = 0",
          "[lancet][cbdg][GraphComplexity][GraphEntanglementIndex]") {
  using lancet::cbdg::GraphComplexity;

  GraphComplexity const cplx;
  CHECK(cplx.GraphEntanglementIndex() == Catch::Approx(0.0).margin(1e-10));
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  GraphMetrics::FormatVcfValue                                            ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("GraphMetrics::FormatVcfValue: 3 comma-separated values",
          "[lancet][cbdg][GraphComplexity]") {
  lancet::caller::GraphMetrics gmtx;
  gmtx.mGraphEntanglementIndex = 2.345;
  gmtx.mTipToPathCovRatio = 0.15;
  gmtx.mMaxSingleDirDegree = 4;

  auto const vcf_val = gmtx.FormatVcfValue();
  // 3 values → 2 commas
  CHECK(std::count(vcf_val.begin(), vcf_val.end(), ',') == 2);
  // Ends with ",4" (the MaxSingleDirDegree)
  CHECK(vcf_val.find(",4") != std::string::npos);
}
