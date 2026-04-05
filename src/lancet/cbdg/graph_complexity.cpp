#include "lancet/cbdg/graph_complexity.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "lancet/cbdg/graph.h"
#include "lancet/cbdg/node.h"

namespace lancet::cbdg {

auto ComputeGraphComplexity(const Graph& graph, const usize component_id) -> GraphComplexity {
  GraphComplexity cx;

  // Local accumulators — V and E are intermediates for CC, not stored on the class.
  usize num_nodes = 0;
  usize num_edges = 0;
  usize unitig_nodes = 0;

  // Accumulators for coverage statistics
  std::vector<f64> coverages;
  std::vector<f64> tip_coverages;
  std::vector<f64> unitig_coverages;

  for (const auto& [nid, node_ptr] : graph.Nodes()) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (node_ptr->GetComponentId() != component_id) continue;
    num_nodes++;

    const auto dflt_sign = node_ptr->SignFor(Kmer::Ordering::DEFAULT);
    usize dflt_dir_edges = 0;
    usize oppo_dir_edges = 0;

    for (const Edge& edge : *node_ptr) {
      if (edge.SrcSign() == dflt_sign) {
        dflt_dir_edges++;
      } else {
        oppo_dir_edges++;
      }
    }

    // Total edges (including mirrors stored at both endpoints); halved below
    num_edges += dflt_dir_edges + oppo_dir_edges;
    const auto max_dir = std::max(dflt_dir_edges, oppo_dir_edges);
    cx.mMaxSingleDirDegree = std::max(cx.mMaxSingleDirDegree, max_dir);

    const bool is_branch = (dflt_dir_edges >= 2 || oppo_dir_edges >= 2);
    if (is_branch) {
      cx.mNumBranchPoints++;
    }

    // Unitig ratio: nodes with exactly 1-in and 1-out (linear chain)
    if (dflt_dir_edges == 1 && oppo_dir_edges == 1) {
      unitig_nodes++;
    }

    // Coverage for CV calculation
    const auto cov = static_cast<f64>(node_ptr->TotalReadSupport());
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
  cx.mCyclomaticComplexity = (num_edges >= num_nodes) ? (num_edges - num_nodes + 1) : 0;

  // Unitig ratio
  cx.mUnitigRatio = num_nodes > 0
      ? static_cast<f64>(unitig_nodes) / static_cast<f64>(num_nodes) : 0.0;

  // Coverage CV (σ/μ)
  if (!coverages.empty()) {
    const f64 sum = std::accumulate(coverages.begin(), coverages.end(), 0.0);
    const f64 mean = sum / static_cast<f64>(coverages.size());
    if (mean > 0.0) {
      f64 sq_sum = 0.0;
      for (const f64 c : coverages) {
        const f64 diff = c - mean;
        sq_sum += diff * diff;
      }
      cx.mCoverageCv = std::sqrt(sq_sum / static_cast<f64>(coverages.size())) / mean;
    }
  }

  // Tip-to-path coverage ratio
  if (!tip_coverages.empty() && !unitig_coverages.empty()) {
    const f64 tip_mean =
        std::accumulate(tip_coverages.begin(), tip_coverages.end(), 0.0) /
        static_cast<f64>(tip_coverages.size());
    const f64 unitig_mean =
        std::accumulate(unitig_coverages.begin(), unitig_coverages.end(), 0.0) /
        static_cast<f64>(unitig_coverages.size());
    cx.mTipToPathCovRatio = unitig_mean > 0.0 ? tip_mean / unitig_mean : 0.0;
  }

  return cx;
}

}  // namespace lancet::cbdg
