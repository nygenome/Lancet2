#include "lancet/cbdg/max_flow.h"

#include "lancet/base/assert.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/graph.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/traversal_index.h"

#include "absl/container/inlined_vector.h"

#include <absl/container/chunked_queue.h>
#include <algorithm>
#include <optional>
#include <string>
#include <vector>

#include <cstddef>

namespace lancet::cbdg {

MaxFlow::MaxFlow(Graph::NodeTable const* graph, NodeIDPair const& /*src_and_snk*/,
                 usize const currk, TraversalIndex const* trav_idx)
    : mGraph(graph), mIndex(trav_idx), mCurrentK(currk) {
  LANCET_ASSERT(mGraph != nullptr)
  LANCET_ASSERT(mIndex != nullptr)
}

// ============================================================================
//  ReconstructWalk — Follow parent pointers from leaf to root
// ============================================================================
//
// The walk tree stores edges as parent-linked arena nodes. To get the full
// ordered edge list (source → ... → sink), we walk backwards from the leaf
// to root (NO_PARENT), collecting edges, then reverse.
//
//   Arena: [0: edge₀,NO_PARENT] [1: edge₁,0] [2: edge₂,1] [3: edge₃,2]
//   Leaf = 3: traverse 3→2→1→0 → edges = [e₃,e₂,e₁,e₀] → reverse → [e₀,e₁,e₂,e₃]
//
auto MaxFlow::ReconstructWalk(std::vector<WalkTreeNode> const& arena, u32 const leaf_idx) const
    -> Walk {
  Walk edges;
  u32 idx = leaf_idx;
  while (idx != TraversalIndex::NO_PARENT) {
    edges.push_back(mIndex->mOrigEdges[arena[idx].mEdgeOrdinal]);
    idx = arena[idx].mParentIdx;
  }
  std::ranges::reverse(edges);
  return edges;
}

// ============================================================================
//  BuildSequence — Assemble haplotype DNA sequence from a walk
// ============================================================================
//
// Each edge in the walk connects two nodes whose sequences overlap by (k-1)
// bases. The haplotype sequence is built by concatenating the non-overlapping
// suffix of each destination node's sequence.
//
auto MaxFlow::BuildSequence(WalkView const walk) const -> Result {
  LANCET_ASSERT(!walk.empty())

  Graph::Path path;
  usize total_seq_len = 0;
  std::vector<std::string> uniq_seqs;
  uniq_seqs.reserve(walk.size() + 1);

  constexpr auto DEFAULT_ORDER = Kmer::Ordering::DEFAULT;
  constexpr auto OPPOSITE_ORDER = Kmer::Ordering::OPPOSITE;
  auto ordering = walk[0].SrcSign() == Kmer::Sign::PLUS ? DEFAULT_ORDER : OPPOSITE_ORDER;

  for (auto const& conn : walk) {
    if (uniq_seqs.empty()) {
      auto const src_itr = mGraph->find(conn.SrcId());
      LANCET_ASSERT(src_itr != mGraph->end())
      LANCET_ASSERT(src_itr->second != nullptr)
      uniq_seqs.emplace_back(src_itr->second->SequenceFor(ordering));
      total_seq_len += uniq_seqs.back().length();
      path.AddNodeCoverage(src_itr->second->TotalReadSupport());
    }

    auto const dst_itr = mGraph->find(conn.DstId());
    LANCET_ASSERT(dst_itr != mGraph->end())
    LANCET_ASSERT(dst_itr->second != nullptr)

    ordering = conn.DstSign() == Kmer::Sign::PLUS ? DEFAULT_ORDER : OPPOSITE_ORDER;
    auto const dst_seq = dst_itr->second->SequenceFor(ordering);
    auto const uniq_seq_len = dst_seq.size() - mCurrentK + 1;
    uniq_seqs.emplace_back(dst_seq.substr(mCurrentK - 1, uniq_seq_len));
    total_seq_len += uniq_seqs.back().length();
    path.AddNodeCoverage(dst_itr->second->TotalReadSupport());
  }

  if (uniq_seqs.empty())
    return std::nullopt;

  path.ReserveSequence(total_seq_len);
  for (auto const& item : uniq_seqs) {
    path.AppendSequence(item);
  }

  LANCET_ASSERT(path.Sequence().length() == total_seq_len)

  path.Finalize();
  return path;
}

// ============================================================================
//  NextPath — Find the Next Walk with At Least One New Edge
// ============================================================================
//
// This is the core walk enumeration routine, matching the original Lancet
// algorithm's semantics but using the walk-tree arena to eliminate the
// exponential walk-vector copying that was the original bottleneck.
//
// ALGORITHM
// ----------
// 1. BFS from source, building a walk tree in an arena.
// 2. Each arena node tracks its accumulated "score" — the count of edges
//    on its walk that are NOT in mTraversedOrdinals.
// 3. When BFS reaches the sink:
//    a. If score > 0 → this walk has at least one new edge. Accept it.
//    b. If score == 0 → this walk only uses already-traversed edges. Skip.
// 4. After accepting a walk, its edges are added to mTraversedOrdinals
//    so subsequent calls will seek walks containing OTHER new edges.
// 5. BFS is bounded by DEFAULT_GRAPH_TRAVERSAL_LIMIT visits to prevent
//    combinatorial blowup on pathological graphs.
//
//   SCORE PROPAGATION IN THE WALK TREE
//   ┌──────────────────────────────────────────────────────┐
//   │ Source outgoing edges: e₀ (new, score=1), e₁ (old)   │
//   │                                                      │
//   │ Arena[0]: e₀, score=1, parent=NO_PARENT              │
//   │ Arena[1]: e₁, score=0, parent=NO_PARENT              │
//   │                                                      │
//   │ Expand Arena[0] → child e₂ (old):                    │
//   │ Arena[2]: e₂, score=1, parent=0   (inherits 1+0)     │
//   │                                                      │
//   │ Expand Arena[1] → child e₃ (new):                    │
//   │ Arena[3]: e₃, score=1, parent=1   (inherits 0+1)     │
//   │                                                      │
//   │ If Arena[2] reaches sink: score=1 > 0 → ACCEPT       │
//   │ Walk = reconstruct(Arena,2) → [e₀, e₂]               │
//   │ Mark e₀ as traversed.                                │
//   └──────────────────────────────────────────────────────┘
//
// WHY SCORE > 0 IS REQUIRED
// ---------------------------
// Without the score check, the algorithm would keep returning the same
// walk (or walks with only already-traversed edges) indefinitely.
// The score ensures monotonic progress: each call returns a walk with
// at least one edge not seen in any previous walk. When no such walk
// exists, the enumeration terminates.
//
auto MaxFlow::NextPath() -> Result {
  std::vector<WalkTreeNode> arena;
  arena.reserve(static_cast<std::size_t>(mIndex->NumNodes()) * 2);

  // BFS frontier: arena indices of tree nodes to expand.
  absl::chunked_queue<u32, 256, 1024> frontier;

  // Seed: outgoing edges from source state.
  EnqueueOutgoingEdges(mIndex->mSrcState, TraversalIndex::NO_PARENT, 0, arena, frontier);

  u32 nvisits = 0;
  std::optional<u32> best_leaf;

  while (!frontier.empty()) {
    nvisits++;
    if (nvisits > Graph::DEFAULT_GRAPH_TRAVERSAL_LIMIT) {
      break;
    }

    u32 const arena_idx = frontier.front();
    frontier.pop_front();
    auto const& node = arena[arena_idx];

    // --- Sink reached: check if this walk has any new edges ---
    if (mIndex->IsSinkState(node.mDstState)) {
      if (node.mScore == 0) {
        continue;  // Only traversed edges → skip this walk
      }

      // Accept first walk with score > 0. Because BFS explores by path length,
      // and we natively rank branches descending by Read Support coverage locally,
      // this guarantees discovering the most biologically prevalent topology first!
      best_leaf = arena_idx;
      break;
    }

    // --- Expand: outgoing edges from this state ---
    EnqueueOutgoingEdges(node.mDstState, arena_idx, node.mScore, arena, frontier);
  }

  // No walk with any new edge found → enumeration complete
  if (!best_leaf.has_value())
    return std::nullopt;

  // Reconstruct the walk and mark its edges as traversed.
  // Walk the arena parent chain to collect ordinals directly — no linear scan.
  Walk const path = ReconstructWalk(arena, *best_leaf);
  u32 mark_idx = *best_leaf;
  while (mark_idx != TraversalIndex::NO_PARENT) {
    mTraversedOrdinals.insert(arena[mark_idx].mEdgeOrdinal);
    mark_idx = arena[mark_idx].mParentIdx;
  }

  return BuildSequence(path);
}

// ============================================================================
//  EnqueueOutgoingEdges — Sort and Push Topologically Dominant Branches
// ============================================================================
//
// Extracts all outgoing edges for a given traversal state, sorts them descending
// by their destination node's read-support coverage, and enqueues them into the
// BFS frontier in a two-pass priority scheme.
//
// PRIORITY:
// 1. Untraversed edges: Increases the walk score, maximizing newly discovered loops.
// 2. Traversed edges: Does not increase score, serves only to connect novel segments.
//
// SORTING HEURISTIC:
// By resolving ties through coverage strength, the BFS naturally explores the most
// biologically dominant allele pathways first. This prevents rare artifacts from
// generating structurally impossible chimeras in the final MSA geometry.
//
void MaxFlow::EnqueueOutgoingEdges(u32 const state_idx, u32 const parent_ai, u32 const parent_score,
                                   std::vector<WalkTreeNode>& arena,
                                   absl::chunked_queue<u32, 256, 1024>& frontier) const {
  auto const& range = mIndex->mAdjRanges[state_idx];
  if (range.mCount == 0) {
    return;
  }

  // Materialize edges into stack array to sort them by destination read support
  absl::InlinedVector<TraversalIndex::OutEdge, 8> out_edges;
  out_edges.reserve(range.mCount);
  for (u32 i = 0; i < range.mCount; i++) {
    out_edges.push_back(mIndex->mAdjList[range.mStart + i]);
  }

  // Sort edges descending by the destination node's TotalReadSupport to organically traverse
  // high-confidence paths first
  std::ranges::sort(
      out_edges,
      [this](TraversalIndex::OutEdge const& lhs, TraversalIndex::OutEdge const& rhs) -> bool {
        auto const lhs_cov =
            this->mIndex->mNodes[TraversalIndex::NodeIdxOf(lhs.mDstState)]->TotalReadSupport();
        auto const rhs_cov =
            this->mIndex->mNodes[TraversalIndex::NodeIdxOf(rhs.mDstState)]->TotalReadSupport();
        return lhs_cov > rhs_cov;
      });

  // Pass 1: untraversed edges (high priority — extend walk score)
  for (auto const& out : out_edges) {
    if (mTraversedOrdinals.contains(out.mEdgeOrdinal)) {
      continue;
    }

    // For source edges, parent_score is technically 0, but since this edge is untraversed, score
    // is 1.
    auto const child_ai = static_cast<u32>(arena.size());
    arena.emplace_back(out.mEdgeOrdinal, out.mDstState, parent_ai, parent_score + 1);
    frontier.push_back(child_ai);
  }

  // Pass 2: already-traversed edges (low priority — no score increase)
  for (auto const& out : out_edges) {
    if (!mTraversedOrdinals.contains(out.mEdgeOrdinal)) {
      continue;
    }

    auto const child_ai = static_cast<u32>(arena.size());
    arena.emplace_back(out.mEdgeOrdinal, out.mDstState, parent_ai, parent_score);
    frontier.push_back(child_ai);
  }
}

}  // namespace lancet::cbdg
