#ifndef SRC_LANCET_CBDG_MAX_FLOW_H_
#define SRC_LANCET_CBDG_MAX_FLOW_H_

#include <deque>
#include <optional>
#include <string>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/container/inlined_vector.h"
#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/graph.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"

namespace lancet::cbdg {

class MaxFlow {
 public:
  explicit MaxFlow(const Graph::NodeTable* graph, const NodeIDPair& src_and_snk, usize currk);

  using Result = std::optional<std::string>;
  [[nodiscard]] auto NextPath() -> Result;

 private:
  static constexpr usize INLINE_EDGES = 8;
  absl::flat_hash_set<Edge> mTraversed;
  absl::InlinedVector<Edge, INLINE_EDGES> mWalkableEdges;

  const Graph::NodeTable* mGraph = nullptr;
  const Node* mSource = nullptr;
  const Node* mSink = nullptr;
  usize mCurrentK = 0;

  using Walk = std::vector<Edge>;
  using WalkView = absl::Span<const Edge>;
  using CandidateWalks = std::deque<Walk>;

  [[nodiscard]] auto BuildNextWalk() -> std::optional<Walk>;

  [[nodiscard]] auto BuildSequence(WalkView walk) const -> Result;
  void PopulateWalkableEdgesInDirection(const Node* src, Kmer::Sign dir);
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_MAX_FLOW_H_
