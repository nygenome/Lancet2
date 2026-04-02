#include "lancet/cbdg/max_flow.h"

#include <algorithm>
#include <optional>
#include <utility>

#include "absl/container/chunked_queue.h"
#include "absl/strings/str_cat.h"
#include "lancet/base/assert.h"

namespace lancet::cbdg {

MaxFlow::MaxFlow(const Graph::NodeTable *graph, const NodeIDPair &src_and_snk,
                 const usize currk, const TraversalIndex *trav_idx)
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
auto MaxFlow::ReconstructWalk(const std::vector<WalkTreeNode> &arena, const u32 leaf_idx) const -> Walk {
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
auto MaxFlow::BuildSequence(const WalkView walk) const -> Result {
  LANCET_ASSERT(!walk.empty())

  usize total_seq_len = 0;
  std::vector<std::string> uniq_seqs;
  uniq_seqs.reserve(walk.size() + 1);

  constexpr auto dflt_order = Kmer::Ordering::DEFAULT;
  constexpr auto oppo_order = Kmer::Ordering::OPPOSITE;
  auto ordering = walk[0].SrcSign() == Kmer::Sign::PLUS ? dflt_order : oppo_order;

  for (const auto &conn : walk) {
    if (uniq_seqs.empty()) {
      const auto src_itr = mGraph->find(conn.SrcId());
      LANCET_ASSERT(src_itr != mGraph->end())
      LANCET_ASSERT(src_itr->second != nullptr)
      uniq_seqs.emplace_back(src_itr->second->SequenceFor(ordering));
      total_seq_len += uniq_seqs.back().length();
    }

    const auto dst_itr = mGraph->find(conn.DstId());
    LANCET_ASSERT(dst_itr != mGraph->end())
    LANCET_ASSERT(dst_itr->second != nullptr)

    ordering = conn.DstSign() == Kmer::Sign::PLUS ? dflt_order : oppo_order;
    const auto dst_seq = dst_itr->second->SequenceFor(ordering);
    const auto uniq_seq_len = dst_seq.size() - mCurrentK + 1;
    uniq_seqs.emplace_back(dst_seq.substr(mCurrentK - 1, uniq_seq_len));
    total_seq_len += uniq_seqs.back().length();
  }

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (uniq_seqs.empty()) return std::nullopt;

  std::string merged_seq;
  merged_seq.reserve(total_seq_len);
  for (const auto &item : uniq_seqs) {
    absl::StrAppend(&merged_seq, item);
  }

  LANCET_ASSERT(merged_seq.length() == total_seq_len)
  return merged_seq;
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
//   │ Source outgoing edges: e₀ (new, score=1), e₁ (old)  │
//   │                                                      │
//   │ Arena[0]: e₀, score=1, parent=NO_PARENT              │
//   │ Arena[1]: e₁, score=0, parent=NO_PARENT              │
//   │                                                      │
//   │ Expand Arena[0] → child e₂ (old):                    │
//   │ Arena[2]: e₂, score=1, parent=0   (inherits 1+0)    │
//   │                                                      │
//   │ Expand Arena[1] → child e₃ (new):                    │
//   │ Arena[3]: e₃, score=1, parent=1   (inherits 0+1)    │
//   │                                                      │
//   │ If Arena[2] reaches sink: score=1 > 0 → ACCEPT      │
//   │ Walk = reconstruct(Arena,2) → [e₀, e₂]              │
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
  arena.reserve(mIndex->NumNodes() * 2);

  // BFS frontier: arena indices of tree nodes to expand.
  absl::chunked_queue<u32, 256, 1024> frontier;

  // Seed: outgoing edges from source state.
  // Enqueue untraversed edges first for BFS priority.
  const auto &src_range = mIndex->mAdjRanges[mIndex->mSrcState];

  // Pass 1: untraversed edges from source (high priority)
  for (u32 i = 0; i < src_range.mCount; i++) {
    const auto &out = mIndex->mAdjList[src_range.mStart + i];
    if (mTraversedOrdinals.contains(out.mEdgeOrdinal)) continue;
    const auto ai = static_cast<u32>(arena.size());
    arena.emplace_back(out.mEdgeOrdinal, out.mDstState, TraversalIndex::NO_PARENT, 1);
    frontier.push_back(ai);
  }

  // Pass 2: already-traversed edges from source (low priority)
  for (u32 i = 0; i < src_range.mCount; i++) {
    const auto &out = mIndex->mAdjList[src_range.mStart + i];
    if (!mTraversedOrdinals.contains(out.mEdgeOrdinal)) continue;
    const auto ai = static_cast<u32>(arena.size());
    arena.emplace_back(out.mEdgeOrdinal, out.mDstState, TraversalIndex::NO_PARENT, 0);
    frontier.push_back(ai);
  }

  u32 nvisits = 0;
  std::optional<u32> best_leaf;
  u32 best_score = 0;

  while (!frontier.empty()) {
    nvisits++;
    if (nvisits > Graph::DEFAULT_GRAPH_TRAVERSAL_LIMIT) break;

    const u32 ai = frontier.front();
    frontier.pop_front();
    const auto &node = arena[ai];

    // --- Sink reached: check if this walk has any new edges ---
    if (mIndex->IsSinkState(node.mDstState)) {
      if (node.mScore > 0 && node.mScore > best_score) {
        best_leaf = ai;
        best_score = node.mScore;
        // Accept first walk with score > 0 (BFS guarantees shortest)
        break;
      }
      // Score == 0: only traversed edges → skip this walk
      continue;
    }

    // --- Expand: outgoing edges from this state ---
    // Enqueue untraversed edges first (two-pass) to match the original code's
    // priority logic without per-node sorting or allocation overhead.
    // No per-walk cycle detection: the 1M visit cap (DEFAULT_GRAPH_TRAVERSAL_LIMIT)
    // bounds exploration, and BFS FIFO order naturally prefers shorter paths.
    const auto &range = mIndex->mAdjRanges[node.mDstState];

    // Pass 1: untraversed edges (high priority — extend walk score)
    for (u32 i = 0; i < range.mCount; i++) {
      const auto &out = mIndex->mAdjList[range.mStart + i];
      if (mTraversedOrdinals.contains(out.mEdgeOrdinal)) continue;
      const auto child_ai = static_cast<u32>(arena.size());
      arena.emplace_back(out.mEdgeOrdinal, out.mDstState, ai, node.mScore + 1);
      frontier.push_back(child_ai);
    }

    // Pass 2: already-traversed edges (low priority — no score increase)
    for (u32 i = 0; i < range.mCount; i++) {
      const auto &out = mIndex->mAdjList[range.mStart + i];
      if (!mTraversedOrdinals.contains(out.mEdgeOrdinal)) continue;
      const auto child_ai = static_cast<u32>(arena.size());
      arena.emplace_back(out.mEdgeOrdinal, out.mDstState, ai, node.mScore);
      frontier.push_back(child_ai);
    }
  }

  // No walk with any new edge found → enumeration complete
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (!best_leaf.has_value()) return std::nullopt;

  // Reconstruct the walk and mark its edges as traversed.
  // Walk the arena parent chain to collect ordinals directly — no linear scan.
  Walk path = ReconstructWalk(arena, *best_leaf);
  u32 mark_idx = *best_leaf;
  while (mark_idx != TraversalIndex::NO_PARENT) {
    mTraversedOrdinals.insert(arena[mark_idx].mEdgeOrdinal);
    mark_idx = arena[mark_idx].mParentIdx;
  }

  return BuildSequence(path);
}

}  // namespace lancet::cbdg
