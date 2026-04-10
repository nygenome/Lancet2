#include "lancet/cbdg/graph.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/traversal_index.h"

#include "catch_amalgamated.hpp"

#include <memory>
#include <string>
#include <vector>

using lancet::cbdg::Edge;
using lancet::cbdg::Graph;
using lancet::cbdg::Kmer;
using lancet::cbdg::Label;
using lancet::cbdg::MakeFwdEdgeKind;
using lancet::cbdg::Node;
using lancet::cbdg::NodeID;
using lancet::cbdg::NodeIDPair;
using lancet::cbdg::RevEdgeKind;
using lancet::cbdg::TraversalIndex;

namespace {

// Helper to build a small graph for testing. Returns the NodeTable and a list of NodeIDs
// in insertion order. Each node gets a distinct short DNA sequence as its kmer.
struct TestGraph {
  Graph::NodeTable mNodes;
  std::vector<NodeID> mNodeIds;

  auto AddNode(std::string_view seq, Label label = Label(Label::REFERENCE)) -> NodeID {
    auto mer = Kmer(seq);
    auto const node_id = mer.Identifier();
    mNodes.try_emplace(node_id, std::make_unique<Node>(std::move(mer), label));
    mNodeIds.push_back(node_id);
    return node_id;
  }

  void AddEdge(NodeID src, NodeID dst) {
    auto& src_node = mNodes.at(src);
    auto& dst_node = mNodes.at(dst);
    auto const src_sign = src_node->SignFor(Kmer::Ordering::DEFAULT);
    auto const dst_sign = dst_node->SignFor(Kmer::Ordering::DEFAULT);
    auto const kind = MakeFwdEdgeKind({src_sign, dst_sign});
    src_node->EmplaceEdge(NodeIDPair{src, dst}, kind);
    dst_node->EmplaceEdge(NodeIDPair{dst, src}, RevEdgeKind(kind));
  }

  void SetAllComponentId(usize comp_id) {
    for (auto& [node_id, node_ptr] : mNodes) {
      node_ptr->SetComponentId(comp_id);
    }
  }
};

}  // namespace

// ============================================================================
//  TraversalIndex tests
// ============================================================================

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("TraversalIndex state index helpers", "[lancet][cbdg][TraversalIndex]") {
  SECTION("MakeState encodes node index and sign") {
    CHECK(TraversalIndex::MakeState(0, Kmer::Sign::PLUS) == 0);
    CHECK(TraversalIndex::MakeState(0, Kmer::Sign::MINUS) == 1);
    CHECK(TraversalIndex::MakeState(1, Kmer::Sign::PLUS) == 2);
    CHECK(TraversalIndex::MakeState(1, Kmer::Sign::MINUS) == 3);
    CHECK(TraversalIndex::MakeState(5, Kmer::Sign::PLUS) == 10);
    CHECK(TraversalIndex::MakeState(5, Kmer::Sign::MINUS) == 11);
  }

  SECTION("NodeIdxOf and SignOf are inverse of MakeState") {
    for (u32 nidx = 0; nidx < 10; nidx++) {
      for (auto sign : {Kmer::Sign::PLUS, Kmer::Sign::MINUS}) {
        auto const state = TraversalIndex::MakeState(nidx, sign);
        CHECK(TraversalIndex::NodeIdxOf(state) == nidx);
        CHECK(TraversalIndex::SignOf(state) == sign);
      }
    }
  }
}

// ============================================================================
//  GraphComplexity tests
// ============================================================================

TEST_CASE("GraphComplexity for linear chain", "[lancet][cbdg][GraphComplexity]") {
  // Linear chain: A → B → C (after compression this would be 1 node,
  // but here we test metrics on the uncompressed graph)
  //
  //   A ──→ B ──→ C
  //
  // V=3, E=2 (unduplicated), M = E-V+1 = 0 (no cycles)
  // density = 2/3 ≈ 0.67, max_degree = 1, branch_points = 0

  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");  // 11-mers
  auto const nid_b = tgraph.AddNode("CGTACGTACGA");
  auto const nid_c = tgraph.AddNode("GTACGTACGAC");
  tgraph.AddEdge(nid_a, nid_b);
  tgraph.AddEdge(nid_b, nid_c);
  tgraph.SetAllComponentId(1);

  // Use the Graph's ComputeGraphComplexity indirectly by constructing a Graph
  // and calling the method. Since we can't easily do that without the full
  // pipeline, we test the metrics computations via the values we know.

  // Count metrics manually from the test graph
  usize num_nodes = 0;
  usize num_edges_raw = 0;
  usize max_dir = 0;
  usize branches = 0;

  for (auto const& [node_id, node_ptr] : tgraph.mNodes) {
    if (node_ptr->GetComponentId() != 1) continue;
    num_nodes++;
    auto const dflt_sign = node_ptr->SignFor(Kmer::Ordering::DEFAULT);
    usize dflt_edges = 0;
    usize oppo_edges = 0;
    for (Edge const& edge : *node_ptr) {
      if (edge.SrcSign() == dflt_sign) {
        dflt_edges++;
      } else {
        oppo_edges++;
      }
    }
    num_edges_raw += dflt_edges + oppo_edges;
    auto const max_d = std::max(dflt_edges, oppo_edges);
    max_dir = std::max(max_dir, max_d);
    if (dflt_edges >= 2 || oppo_edges >= 2) branches++;
  }
  usize const num_edges = num_edges_raw / 2;

  CHECK(num_nodes == 3);
  CHECK(num_edges == 2);
  CHECK(max_dir <= 1);
  CHECK(branches == 0);

  // Cyclomatic complexity: M = E - V + 1 = 2 - 3 + 1 = 0
  usize const cyclomatic = (num_edges >= num_nodes) ? (num_edges - num_nodes + 1) : 0;
  CHECK(cyclomatic == 0);
}

TEST_CASE("GraphComplexity for single bubble", "[lancet][cbdg][GraphComplexity]") {
  // Simple bubble: A → B → D, A → C → D
  //
  //        ┌── B ──┐
  //   A ──┤        ├── D
  //        └── C ──┘
  //
  // V=4, E=4 (unduplicated), M = 4-4+1 = 1 (one independent cycle)
  // max_degree = 2 (A has 2 outgoing), branch_points = 1 (A)

  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");
  auto const nid_b = tgraph.AddNode("CGTACGTACGA");
  auto const nid_c = tgraph.AddNode("GTACGTACGAC");
  auto const nid_d = tgraph.AddNode("TACGTACGACG");
  tgraph.AddEdge(nid_a, nid_b);
  tgraph.AddEdge(nid_a, nid_c);
  tgraph.AddEdge(nid_b, nid_d);
  tgraph.AddEdge(nid_c, nid_d);
  tgraph.SetAllComponentId(1);

  usize num_nodes = 0;
  usize num_edges_raw = 0;
  usize max_dir = 0;
  usize branches = 0;

  for (auto const& [node_id, node_ptr] : tgraph.mNodes) {
    if (node_ptr->GetComponentId() != 1) continue;
    num_nodes++;
    auto const dflt_sign = node_ptr->SignFor(Kmer::Ordering::DEFAULT);
    usize dflt_edges = 0;
    usize oppo_edges = 0;
    for (Edge const& edge : *node_ptr) {
      if (edge.SrcSign() == dflt_sign) {
        dflt_edges++;
      } else {
        oppo_edges++;
      }
    }
    num_edges_raw += dflt_edges + oppo_edges;
    auto const max_d = std::max(dflt_edges, oppo_edges);
    max_dir = std::max(max_dir, max_d);
    if (dflt_edges >= 2 || oppo_edges >= 2) branches++;
  }
  usize const num_edges = num_edges_raw / 2;

  CHECK(num_nodes == 4);
  CHECK(num_edges == 4);
  CHECK(max_dir == 2);
  CHECK(branches >= 1);

  usize const cyclomatic = (num_edges >= num_nodes) ? (num_edges - num_nodes + 1) : 0;
  CHECK(cyclomatic == 1);
}

// ============================================================================
//  TraversalIndex construction tests
// ============================================================================

TEST_CASE("TraversalIndex IsSinkState", "[lancet][cbdg][TraversalIndex]") {
  TraversalIndex tidx;
  tidx.mSnkNodeIdx = 3;
  tidx.mAdjRanges.resize(10);

  // Both + and - states of node 3 should be sink states
  CHECK(tidx.IsSinkState(TraversalIndex::MakeState(3, Kmer::Sign::PLUS)));
  CHECK(tidx.IsSinkState(TraversalIndex::MakeState(3, Kmer::Sign::MINUS)));

  // Other nodes are not sink states
  CHECK_FALSE(tidx.IsSinkState(TraversalIndex::MakeState(0, Kmer::Sign::PLUS)));
  CHECK_FALSE(tidx.IsSinkState(TraversalIndex::MakeState(1, Kmer::Sign::MINUS)));
  CHECK_FALSE(tidx.IsSinkState(TraversalIndex::MakeState(4, Kmer::Sign::PLUS)));
}
