#ifndef SRC_LANCET_CBDG_TRAVERSAL_INDEX_H_
#define SRC_LANCET_CBDG_TRAVERSAL_INDEX_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"

#include <limits>
#include <vector>

namespace lancet::cbdg {

// ============================================================================
//  TraversalIndex — Flat, Cache-Friendly Adjacency List for Graph Traversal
// ============================================================================
//
// MOTIVATION
// ----------
// The graph's hash-map-based NodeTable (flat_hash_map<NodeID, NodePtr>) is
// optimized for dynamic mutation (insertion, deletion, lookup by hash key).
// However, traversal algorithms (cycle detection, max-flow path finding)
// do not mutate the graph — they only read topology and track per-node state
// (visited/color). Hash-based state tracking is slow in tight loops because
// it causes unpredictable cache misses and expensive hash computations.
//
// APPROACH
// --------
// After all graph mutations are complete (pruning, compression, tip removal),
// we build a read-only TraversalIndex that maps the graph into contiguous
// integer-indexed arrays. All traversal state (DFS colors, BFS visited,
// edge flow flags) becomes simple array operations — O(1) with cache-line-
// friendly sequential access.
//
// NODE-STATE MODEL (Bidirected Graph)
// ------------------------------------
// In the BCALM2 bidirected model, each node has two traversal states:
// one for each sign direction (+ and -). A valid walk must satisfy sign
// continuity: edge_i.DstSign == edge_{i+1}.SrcSign. This means the same
// node reached via '+' and via '-' are different traversal states.
//
// We encode this as a flat state index:
//
//   state_idx = node_flat_idx * 2 + sign_offset
//   sign_offset = 0 for PLUS, 1 for MINUS
//
//   Example for 3 nodes:
//   ┌──────────┬───────────────────────────────┐
//   │ node_idx │  state=0(+)     state=1(-)    │
//   ├──────────┼───────────────────────────────┤
//   │    0     │     0              1          │
//   │    1     │     2              3          │
//   │    2     │     4              5          │
//   └──────────┴───────────────────────────────┘
//
// FLAT ADJACENCY LIST (CSR-like)
// --------------------------------
// For each state, we store a contiguous range of outgoing edges in mAdjList:
//
//   mAdjRanges[state] = { start, count }
//   mAdjList[start .. start+count) = outgoing edges from that state
//
//   ┌──────────┐   ┌────────────────────────────────────────┐
//   │AdjRanges │   │            AdjList (packed)            │
//   │[state_0] │──>│ edge_0  edge_1  │  edge_2  │  edge_3   │
//   │[state_1] │───────────────────────>│        │          │
//   │[state_2] │────────────────────────────────>│          │
//   └──────────┘   └────────────────────────────────────────┘
//
// LIFECYCLE
// ---------
// Built once per (component, k-value) after graph pruning completes.
// Not maintained during graph mutations — discarded and rebuilt if k changes.
// Consumed by HasCycle (3-color DFS) and MaxFlow (Edmonds-Karp BFS).
//
// ============================================================================

class TraversalIndex {
 public:
  /// An outgoing edge in the flat adjacency list.
  struct OutEdge {
    u32 mDstState;     ///< Destination state index (node_idx*2 + sign_offset)
    u32 mEdgeOrdinal;  ///< Index into mOrigEdges for edge identity/reconstruction
  };

  /// Range into mAdjList for one state. [mStart, mStart + mCount)
  struct AdjRange {
    u32 mStart = 0;
    u32 mCount = 0;
  };

  // --- Core flat adjacency data ---
  std::vector<AdjRange> mAdjRanges;  ///< Indexed by state_idx (size = 2 * num_nodes)
  std::vector<OutEdge> mAdjList;     ///< All outgoing edges, packed contiguously by state

  // --- Edge identity for reconstruction and flow tracking ---
  std::vector<Edge> mOrigEdges;  ///< Indexed by edge ordinal. Used to reconstruct walks.

  // --- Node mapping for sequence reconstruction ---
  std::vector<Node const*> mNodes;  ///< Indexed by node_flat_idx (size = V)
  std::vector<NodeID> mNodeIds;     ///< Indexed by node_flat_idx (size = V)

  // --- Source and sink identifiers ---
  u32 mSrcState = 0;    ///< State index of source node in its default sign
  u32 mSnkNodeIdx = 0;  ///< Flat node index of sink (reachable via either sign)

  /// Sentinel value for "no parent" in walk tree nodes.
  static constexpr u32 NO_PARENT = std::numeric_limits<u32>::max();

  // --- Accessors ---
  [[nodiscard]] auto NumNodes() const -> u32 { return static_cast<u32>(mNodes.size()); }
  [[nodiscard]] auto NumStates() const -> u32 { return static_cast<u32>(mAdjRanges.size()); }
  [[nodiscard]] auto NumEdgeOrdinals() const -> u32 { return static_cast<u32>(mOrigEdges.size()); }

  /// Check if a state index corresponds to the sink node (either + or - sign).
  [[nodiscard]] auto IsSinkState(u32 const state_idx) const -> bool {
    return NodeIdxOf(state_idx) == mSnkNodeIdx;
  }

  // --- Static state index helpers ---

  [[nodiscard]] static constexpr auto MakeState(u32 const node_idx, Kmer::Sign const sign) -> u32 {
    return (node_idx * 2) + (sign == Kmer::Sign::PLUS ? 0 : 1);
  }

  [[nodiscard]] static constexpr auto NodeIdxOf(u32 const state_idx) -> u32 {
    return state_idx / 2;
  }

  [[nodiscard]] static constexpr auto SignOf(u32 const state_idx) -> Kmer::Sign {
    return (state_idx % 2 == 0) ? Kmer::Sign::PLUS : Kmer::Sign::MINUS;
  }
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_TRAVERSAL_INDEX_H_
