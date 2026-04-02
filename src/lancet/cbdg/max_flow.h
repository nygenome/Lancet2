#ifndef SRC_LANCET_CBDG_MAX_FLOW_H_
#define SRC_LANCET_CBDG_MAX_FLOW_H_

#include <optional>
#include <string>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/graph.h"

namespace lancet::cbdg {

// ============================================================================
//  MaxFlow — Sequential Walk Enumeration for Haplotype Assembly
// ============================================================================
//
// ALGORITHM OVERVIEW
// -------------------
// Finds all distinct source-to-sink walks in the bidirected de Bruijn graph
// by iteratively calling NextPath(). Each call performs a BFS from source
// that scores candidate walks by the number of NOT-YET-TRAVERSED edges they
// contain. When the BFS finds a walk to sink with score > 0 (at least one
// new edge), it is returned and its edges are marked as traversed. The next
// call will preferentially explore walks using new edges.
//
// This is the ORIGINAL Lancet algorithm, reimplemented with the walk-tree
// arena to eliminate the exponential walk-vector copying and per-walk
// hashing that caused the original performance bottleneck.
//
//   SEQUENTIAL WALK ENUMERATION
//   ┌──────────────────────────────────────────────────────┐
//   │ Call 1: BFS finds walk with highest # of new edges   │
//   │         Return Src→A→M→Sink  (3 new edges)          │
//   │         Mark {Src→A, A→M, M→Sink} as traversed      │
//   │                                                      │
//   │ Call 2: BFS explores all walks. Walks reusing only   │
//   │         traversed edges get score=0 and are skipped. │
//   │         Return Src→B→M→Sink  (1 new edge: Src→B)    │
//   │         Mark {Src→B} as traversed                    │
//   │                                                      │
//   │ Call 3: No walk has any new edge → return nullopt    │
//   │         → enumeration complete                       │
//   └──────────────────────────────────────────────────────┘
//
// WALK TREE ARENA (perf improvement over original)
// -------------------------------------------------
// Instead of BFS that copies entire walk vectors at each extension level
// (the bottleneck in the original code), each BFS node in the arena stores
// only an edge ordinal (4B) + parent index (4B) + accumulated score (4B).
// Full walks are reconstructed by parent-pointer traversal only when a
// walk reaches the sink with score > 0.
//
//   OLD: Walk copy per BFS extension → O(B^L × L) memory
//   NEW: Arena node per BFS extension → O(B^L × 16B) memory, no copies
//
// Both approaches are bounded by DEFAULT_GRAPH_TRAVERSAL_LIMIT (1M BFS
// visits per call) to prevent combinatorial blowup on pathological graphs.
//
// ============================================================================
class MaxFlow {
 public:
  explicit MaxFlow(const Graph::NodeTable* graph, const NodeIDPair& src_and_snk, usize currk, const TraversalIndex* trav_idx);

  using Result = std::optional<std::string>;

  /// Find the next walk from source to sink that contains at least one
  /// edge not yet traversed by any previous walk. Returns nullopt when
  /// no such walk exists (all edges covered or traversal limit reached).
  [[nodiscard]] auto NextPath() -> Result;

 private:
  const Graph::NodeTable* mGraph = nullptr;
  const TraversalIndex* mIndex = nullptr;
  usize mCurrentK = 0;

  /// Set of edge ordinals already traversed by previous walks.
  /// Edges in this set get score 0; walks must have at least one
  /// edge NOT in this set to be accepted.
  absl::flat_hash_set<u32> mTraversedOrdinals;

  using Walk = std::vector<Edge>;
  using WalkView = absl::Span<const Edge>;

  /// Walk tree node for BFS. Stores:
  ///   - mEdgeOrdinal: index into TraversalIndex::mOrigEdges
  ///   - mDstState: state reached (for cycle detection)
  ///   - mParentIdx: back-link in arena (NO_PARENT for root)
  ///   - mScore: accumulated count of un-traversed edges on this walk
  struct WalkTreeNode {
    u32 mEdgeOrdinal;
    u32 mDstState;
    u32 mParentIdx;
    u32 mScore;
  };

  /// Reconstruct the full Walk (edge vector) from a tree leaf back to root.
  [[nodiscard]] auto ReconstructWalk(const std::vector<WalkTreeNode>& arena, u32 leaf_idx) const -> Walk;

  /// Build haplotype sequence string from a completed walk.
  [[nodiscard]] auto BuildSequence(WalkView walk) const -> Result;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_MAX_FLOW_H_
