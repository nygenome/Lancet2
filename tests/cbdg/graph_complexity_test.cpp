#include "lancet/cbdg/graph_complexity.h"

#include <cmath>
#include <vector>

#include "absl/types/span.h"
#include "catch_amalgamated.hpp"
#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  GraphComplexity::IsComplex Thresholds                                  ║
// ╚══════════════════════════════════════════════════════════════════════════╝

// NOTE: We cannot directly construct GraphComplexity with specific field values
// because the fields are private (populated only by ComputeGraphComplexity).
// IsComplex() tests use the default-constructed state plus the public constants.

TEST_CASE("GraphComplexity::IsComplex default is not complex", "[lancet][cbdg][graph_complexity]") {
  using lancet::cbdg::GraphComplexity;

  const GraphComplexity cx;
  CHECK_FALSE(cx.IsComplex());
  CHECK(cx.CyclomaticComplexity() == 0);
  CHECK(cx.NumBranchPoints() == 0);
}

TEST_CASE("GraphComplexity::IsComplex thresholds are correct", "[lancet][cbdg][graph_complexity]") {
  using lancet::cbdg::GraphComplexity;

  CHECK(GraphComplexity::MAX_CYCLOMATIC_COMPLEXITY == 50);
  CHECK(GraphComplexity::MAX_BRANCH_POINTS == 50);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  GraphComplexity::GraphEntanglementIndex (GEI)                          ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("GEI: default-constructed → GEI = 0", "[lancet][cbdg][graph_complexity][gei]") {
  using lancet::cbdg::GraphComplexity;

  const GraphComplexity cx;
  CHECK(cx.GraphEntanglementIndex() == Catch::Approx(0.0).margin(1e-10));
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  GraphMetrics::FormatVcfValue                                           ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("GraphMetrics::FormatVcfValue: 3 comma-separated values", "[lancet][cbdg][graph_complexity][format]") {
  lancet::caller::RawVariant::GraphMetrics gm;
  gm.graph_entanglement_index = 2.345;
  gm.tip_to_path_cov_ratio = 0.15;
  gm.max_single_dir_degree = 4;

  const auto vcf = gm.FormatVcfValue();
  // 3 values → 2 commas
  CHECK(std::count(vcf.begin(), vcf.end(), ',') == 2);
  // Ends with ",4" (the MaxSingleDirDegree)
  CHECK(vcf.find(",4") != std::string::npos);
}
