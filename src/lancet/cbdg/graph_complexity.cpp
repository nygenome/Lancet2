#include "lancet/cbdg/graph_complexity.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/graph.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include <cmath>

namespace lancet::cbdg {

auto ComputeGraphComplexity(Graph const& graph, usize const component_id) -> GraphComplexity {
  GraphComplexity cplx;

  // Local accumulators — V and E are intermediates for CC, not stored on the class.
  usize num_nodes = 0;
  usize num_edges = 0;
  usize unitig_nodes = 0;

  // Accumulators for coverage statistics
  std::vector<f64> coverages;
  std::vector<f64> tip_coverages;
  std::vector<f64> unitig_coverages;

  for (auto const& [nid, node_ptr] : graph.Nodes()) {
    if (node_ptr->GetComponentId() != component_id)
      continue;
    num_nodes++;

    auto const dflt_sign = node_ptr->SignFor(Kmer::Ordering::DEFAULT);
    usize dflt_dir_edges = 0;
    usize oppo_dir_edges = 0;

    for (Edge const& edge : *node_ptr) {
      if (edge.SrcSign() == dflt_sign) {
        dflt_dir_edges++;
      } else {
        oppo_dir_edges++;
      }
    }

    // Total edges (including mirrors stored at both endpoints); halved below
    num_edges += dflt_dir_edges + oppo_dir_edges;
    auto const max_dir = std::max(dflt_dir_edges, oppo_dir_edges);
    cplx.mMaxSingleDirDegree = std::max(cplx.mMaxSingleDirDegree, max_dir);

    bool const is_branch = (dflt_dir_edges >= 2 || oppo_dir_edges >= 2);
    if (is_branch) {
      cplx.mNumBranchPoints++;
    }

    // Unitig ratio: nodes with exactly 1-in and 1-out (linear chain)
    if (dflt_dir_edges == 1 && oppo_dir_edges == 1) {
      unitig_nodes++;
    }

    // Coverage for CV calculation
    auto const cov = static_cast<f64>(node_ptr->TotalReadSupport());
    coverages.push_back(cov);

    // Categorize for tip-to-path ratio
    if (dflt_dir_edges == 0 || oppo_dir_edges == 0) {
      tip_coverages.push_back(cov);
    } else if (dflt_dir_edges == 1 && oppo_dir_edges == 1) {
      unitig_coverages.push_back(cov);
    }
  }

  // Each edge stored at both endpoints (forward + mirror) → divide by 2
  num_edges /= 2;
  cplx.mCyclomaticComplexity = (num_edges >= num_nodes) ? (num_edges - num_nodes + 1) : 0;

  // Unitig ratio
  cplx.mUnitigRatio =
      num_nodes > 0 ? static_cast<f64>(unitig_nodes) / static_cast<f64>(num_nodes) : 0.0;

  // Coverage CV (σ/μ)
  if (!coverages.empty()) {
    f64 const sum = std::accumulate(coverages.begin(), coverages.end(), 0.0);
    f64 const mean = sum / static_cast<f64>(coverages.size());
    if (mean > 0.0) {
      f64 sq_sum = 0.0;
      for (f64 const cov_val : coverages) {
        f64 const diff = cov_val - mean;
        sq_sum += diff * diff;
      }
      cplx.mCoverageCv = std::sqrt(sq_sum / static_cast<f64>(coverages.size())) / mean;
    }
  }

  // Tip-to-path coverage ratio
  if (!tip_coverages.empty() && !unitig_coverages.empty()) {
    f64 const tip_mean = std::accumulate(tip_coverages.begin(), tip_coverages.end(), 0.0) /
                         static_cast<f64>(tip_coverages.size());
    f64 const unitig_mean = std::accumulate(unitig_coverages.begin(), unitig_coverages.end(), 0.0) /
                            static_cast<f64>(unitig_coverages.size());
    cplx.mTipToPathCovRatio = unitig_mean > 0.0 ? tip_mean / unitig_mean : 0.0;
  }

  return cplx;
}

}  // namespace lancet::cbdg
