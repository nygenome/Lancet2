#include "lancet/cbdg/graph.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/chunked_queue.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "lancet/base/assert.h"
#include "lancet/base/logging.h"
#include "lancet/base/repeat.h"
#include "lancet/base/sliding.h"
#include "lancet/base/timer.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/max_flow.h"
#include "lancet/cbdg/node.h"
#include "lancet/hts/phred_quality.h"
#include "spdlog/fmt/bundled/core.h"
#include "spdlog/fmt/bundled/ostream.h"

namespace lancet::cbdg {

/// Pipeline architecture for haplotype assembly from a colored de Bruijn graph:
///
///  ┌─────────────┐
///  │ Outer loop:  │  Iterate k from min_k to max_k in steps of mKmerStepLen.
///  │ k-value scan │  If haplotypes are found at any k, stop.
///  │              │  If a cycle is detected or graph is too complex, abandon
///  │              │  this k and continue to next (via should_retry_kmer flag).
///  └──────┬──────┘
///         │
///         ▼
///  ┌─────────────┐
///  │  BuildGraph  │  Build the bidirected de Bruijn graph from reads + reference.
///  │  + prune     │  Remove low-coverage nodes, mark connected components.
///  └──────┬──────┘
///         │
///         ▼
///  ┌─────────────┐   For each connected component with valid source/sink anchors:
///  │ Inner loop:  │
///  │ per-component│   1. Compress linear chains, remove low-cov nodes, remove tips
///  │              │   2. Build TraversalIndex (flat adjacency list) on frozen graph
///  │              │   3. HasCycle via O(V+E) three-color DFS on flat arrays
///  │              │   4. Log graph complexity metrics (cyclomatic, density, etc.)
///  │              │   5. Enumerate all source→sink walks via BFS walk-tree
///  └──────┬──────┘
///         │
///         ▼
///  ┌─────────────┐
///  │  Assemble    │  Deduplicate haplotype sequences, prepend reference anchor.
///  │  + return    │  Return assembled haplotypes for genotyping / variant calling.
///  └─────────────┘
///
/// https://github.com/GATB/bcalm/blob/v2.2.3/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
auto Graph::BuildComponentHaplotypes(RegionPtr region, ReadList reads) -> Result {
  mReads = reads;
  mRegion = std::move(region);

  Timer timer;
  GraphHaps per_comp_haplotypes;
  std::string_view ref_anchor_seq;
  std::vector<usize> anchor_start_idxs;
  absl::flat_hash_set<MateMer> mate_mers;

  static constexpr usize DEFAULT_EST_NUM_NODES = 32768;
  static constexpr usize DEFAULT_MIN_ANCHOR_LENGTH = 150;
  static constexpr f64 DEFAULT_PCT_NODES_NEEDED = 10.0;

  // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
  const auto reg_str = mRegion->ToSamtoolsRegion();
  mCurrK = mParams.mMinKmerLen - mParams.mKmerStepLen;

  // Outer loop: increment k and retry until haplotypes are found or k is exhausted.
  // Replaces the old `goto IncrementKmerAndRetry` with structured control flow.
  while (per_comp_haplotypes.empty() && (mCurrK + mParams.mKmerStepLen) <= mParams.mMaxKmerLen) {
    mCurrK += mParams.mKmerStepLen;
    timer.Reset();
    mSourceAndSinkIds = {0, 0};
    mNodes.reserve(DEFAULT_EST_NUM_NODES);

    // Skip this k if the reference itself has a repeated k-mer — the de Bruijn
    // graph would contain a cycle by construction, making assembly pointless.
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (HasExactOrApproxRepeat(mRegion->SeqView(), mCurrK)) continue;

    mNodes.clear();
    BuildGraph(mate_mers);
    LOG_TRACE("Done building graph for {} with k={}, nodes={}, reads={}", reg_str, mCurrK, mNodes.size(), mReads.size())

    RemoveLowCovNodes(0);
    mNodes.rehash(0);
    WriteDotDevelop(FIRST_LOW_COV_REMOVAL, 0);

    const auto components = MarkConnectedComponents();
    per_comp_haplotypes.reserve(components.size());
    anchor_start_idxs.reserve(components.size());
    LOG_TRACE("Found {} connected components in graph for {} with k={}", components.size(), reg_str, mCurrK)

    // Inner loop: process each connected component with valid source/sink anchors.
    // The should_retry_kmer flag is set to true by cycle detection to abandon all
    // remaining components at this k and retry at a higher k value.
    bool should_retry_kmer = false;
    for (const auto& cinfo : components) {
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (should_retry_kmer) break;
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (cinfo.mPctNodes < DEFAULT_PCT_NODES_NEEDED) continue;

      const auto comp_id = cinfo.mCompId;
      const auto source = FindSource(comp_id);
      const auto sink = FindSink(comp_id);

      if (!source.mFoundAnchor || !sink.mFoundAnchor || source.mAnchorId == sink.mAnchorId) {
        LOG_TRACE("Skipping comp{} in graph for {} because source/sink was not found", comp_id, reg_str)
        continue;
      }

      const auto current_anchor_length = RefAnchorLength(source, sink, mCurrK);
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (current_anchor_length < DEFAULT_MIN_ANCHOR_LENGTH) continue;

      LOG_TRACE("Found {}bp ref anchor for {} comp={} with k={}", current_anchor_length, reg_str, comp_id, mCurrK)

      std::vector<std::string> haplotypes;
      mSourceAndSinkIds = NodeIDPair{source.mAnchorId, sink.mAnchorId};
      ref_anchor_seq = mRegion->SeqView().substr(source.mRefOffset, current_anchor_length);
      WriteDotDevelop(FOUND_REF_ANCHORS, comp_id);

      // Pruning pipeline: compress linear chains, remove low-coverage nodes, remove tips.
      // NOTE: The pre-compression HasCycle call (old line 103) has been removed. Graph
      // compression only merges degree-2 linear chain nodes — it cannot introduce or
      // remove cycles. Therefore cycle detection before compression was redundant with
      // the cycle detection after compression. We now check only once, on the smaller
      // (compressed) graph, which is faster.
      CompressGraph(comp_id);
      WriteDotDevelop(FIRST_COMPRESSION, comp_id);
      RemoveLowCovNodes(comp_id);
      WriteDotDevelop(SECOND_LOW_COV_REMOVAL, comp_id);
      CompressGraph(comp_id);
      WriteDotDevelop(SECOND_COMPRESSION, comp_id);
      RemoveTips(comp_id);
      WriteDotDevelop(SHORT_TIP_REMOVAL, comp_id);

      // Build the flat traversal index on the frozen (fully-pruned) graph.
      // This maps NodeID -> contiguous u32 and constructs the CSR adjacency list.
      // Both HasCycle and MaxFlow operate on this flat structure for O(1) state tracking.
      const auto trav_idx = BuildTraversalIndex(comp_id);

      // O(V+E) cycle detection using three-color DFS on the flat adjacency list.
      // See HasCycle() implementation for bidirected sign-continuity handling.
      if (HasCycle(trav_idx)) {
        LOG_TRACE("Cycle found in pruned graph for {} comp={} with k={}", reg_str, comp_id, mCurrK)
        should_retry_kmer = true;
        break;
      }

      // Log graph complexity metrics for debugging / correlating with runtime.
      // All metrics are O(V+E) to compute and help identify pathological windows.
      // Skip walk enumeration on pathological graphs — retry with larger k to
      // collapse branches. Same control flow as the HasCycle guard above.
      const auto cx = ComputeGraphComplexity(comp_id);
      if (cx.IsComplex()) {
        LOG_DEBUG("Graph too complex for {} comp={} k={}: CC={} branches={}",
                  reg_str, comp_id, mCurrK, cx.mCyclomaticComplexity, cx.mNumBranchPoints)
        should_retry_kmer = true;
        break;
      }

      LOG_DEBUG("Graph complexity stats for {} comp={} k={}: V={} E={} cyclomatic={} branches={}",
                reg_str, comp_id, mCurrK, cx.mNumNodes, cx.mNumEdges, cx.mCyclomaticComplexity, cx.mNumBranchPoints)
      WriteDot(State::FULLY_PRUNED_GRAPH, comp_id);
      LOG_TRACE("Starting walk enumeration for {} with k={}, num_nodes={}", reg_str, mCurrK, mNodes.size())
      MaxFlow max_flow(&mNodes, mSourceAndSinkIds, mCurrK, &trav_idx);
      auto path_seq = max_flow.NextPath();

      while (path_seq) {
        LOG_DEBUG("Assembled {}bp path sequence for {} comp={} with k={}", path_seq->length(), reg_str, comp_id, mCurrK)
        haplotypes.emplace_back(std::move(*path_seq));
        path_seq = max_flow.NextPath();
      }

      if (!haplotypes.empty()) {
        std::ranges::sort(haplotypes);
        const auto dup_range = std::ranges::unique(haplotypes);
        haplotypes.erase(dup_range.begin(), dup_range.end());
        haplotypes.emplace(haplotypes.begin(), ref_anchor_seq);
        per_comp_haplotypes.emplace_back(std::move(haplotypes));
        anchor_start_idxs.emplace_back(source.mRefOffset);
      }
    }

    // If any component triggered a retry, discard partial results and try next k
    if (should_retry_kmer) {
      per_comp_haplotypes.clear();
      anchor_start_idxs.clear();
      continue;
    }
  }

  static const auto summer = [](const u64 sum, const auto& comp_haps) -> u64 { return sum + comp_haps.size() - 1; };
  const auto num_asm_haps = std::accumulate(per_comp_haplotypes.cbegin(), per_comp_haplotypes.cend(), 0, summer);

  // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
  const auto human_rt = timer.HumanRuntime();

  LOG_TRACE("Assembled {} haplotypes for {} with k={} in {}", num_asm_haps, reg_str, mCurrK, human_rt)
  return {.mGraphHaplotypes = per_comp_haplotypes, .mAnchorStartIdxs = anchor_start_idxs};
}

void Graph::CompressGraph(const usize component_id) {
  absl::flat_hash_set<NodeID> remove_nids;
  remove_nids.reserve(mNodes.size());

  for (NodeTable::const_reference item : mNodes) {
    if (item.second->GetComponentId() != component_id) continue;  // NOLINT(readability-braces-around-statements)
    if (remove_nids.contains(item.first)) continue;               // NOLINT(readability-braces-around-statements)

    CompressNode(item.first, Kmer::Ordering::DEFAULT, remove_nids);
    CompressNode(item.first, Kmer::Ordering::OPPOSITE, remove_nids);
  }

  if (!remove_nids.empty()) {
    // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
    const auto region_str = mRegion->ToSamtoolsRegion();
    LOG_TRACE("Compressed {} nodes for {} in comp{} with k={}", remove_nids.size(), region_str, component_id, mCurrK)
    // NOLINTNEXTLINE(readability-braces-around-statements)
    for (const auto nid : remove_nids) RemoveNode(mNodes.find(nid));
  }
}

void Graph::CompressNode(NodeID nid, const Kmer::Ordering ord, NodeIdSet& compressed_ids) const {
  const auto node_itr = mNodes.find(nid);
  LANCET_ASSERT(node_itr != mNodes.end())
  LANCET_ASSERT(node_itr->second != nullptr)

  auto compressible_edge = FindCompressibleEdge(*node_itr->second, ord);
  while (compressible_edge.has_value()) {
    const Edge src2obdy = compressible_edge.value();
    LANCET_ASSERT(src2obdy.SrcId() == nid)
    const auto obdy_itr = mNodes.find(src2obdy.DstId());
    LANCET_ASSERT(obdy_itr != mNodes.end())
    LANCET_ASSERT(obdy_itr->second != nullptr)

    node_itr->second->Merge(*(obdy_itr->second), src2obdy.Kind(), mCurrK);
    node_itr->second->EraseEdge(src2obdy);  // src -->X--> old_buddy

    const auto rev_src2obdy_src_sign = Kmer::RevSign(src2obdy.SrcSign());
    for (const Edge& obdy2nbdy : *(obdy_itr->second)) {
      // Skip if this is old_buddy --> src edge before merging edges
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (obdy2nbdy == src2obdy.MirrorEdge()) continue;

      LANCET_ASSERT(!obdy2nbdy.IsSelfLoop())
      LANCET_ASSERT(obdy2nbdy.DstId() != node_itr->second->Identifier())

      const auto nbdy_itr = mNodes.find(obdy2nbdy.DstId());
      LANCET_ASSERT(nbdy_itr != mNodes.end())
      LANCET_ASSERT(nbdy_itr->second != nullptr)

      // src --> old_buddy --> new_buddy
      // Create src --> new_buddy edge from src --> old_buddy and old_buddy --> new_buddy edges.
      const auto ne_src_sign = src2obdy.DstSign() != obdy2nbdy.SrcSign() ? rev_src2obdy_src_sign : src2obdy.SrcSign();
      const auto src2nbdy = Edge({nid, obdy2nbdy.DstId()}, MakeFwdEdgeKind({ne_src_sign, obdy2nbdy.DstSign()}));

      node_itr->second->EmplaceEdge(src2nbdy);               // src --> new_buddy
      nbdy_itr->second->EmplaceEdge(src2nbdy.MirrorEdge());  // new_buddy --> src
      nbdy_itr->second->EraseEdge(obdy2nbdy.MirrorEdge());   // new_buddy -->X--> old_buddy
    }

    compressed_ids.insert(src2obdy.DstId());
    compressible_edge = FindCompressibleEdge(*node_itr->second, ord);
  }
}

auto Graph::FindCompressibleEdge(const Node& src, const Kmer::Ordering ord) const -> std::optional<Edge> {
  // abc_nbour --> abc_bdy --> src --> xyz_bdy --> xyz_nbour
  // abc_nbour <-- abc_bdy <-- src <-- xyz_bdy --> xyz_nbour
  //
  // Pre-requisites:
  // * Assume we are trying to merge contents of abc_bdy into src. Then result will be src --> abc_bdy edge.
  // * If ord == Default, then result src --> abc_bdy edge must have SrcSign() same as src's default sign
  // * If ord == Opposite, then result src --> abc_bdy edge must have SrcSign() same as src's opposite sign
  // * This expected SrcSign for the result src --> abc_bdy edge is named as `merge_sign`
  //
  // In order for src to be compressible with abc_bdy node, we need to fulfill the following conditions:
  // 1. src must have atmost 2 and atleast 1 outgoing edges and not contain any self loops.
  // 2. src must only have only one out edge where SrcSign is same as merge_sign. i.e. the src --> abc_bdy edge.
  // 3. Another src out edge if present must be in opposite direction of src --> abc_bdy to xyz_bdy.
  // 4. Both src --> abc_bdy && src --> xyz_bdy must be buddy edges from src node.
  // If all these conditions are satisfied, then src --> abc_bdy edge is returned. std::nullopt is returned otherwise.

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (src.NumOutEdges() > 2 || src.NumOutEdges() == 0 || src.HasSelfLoop()) return std::nullopt;

  const auto mergeable_edges = src.FindEdgesInDirection(ord);
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mergeable_edges.size() != 1) return std::nullopt;

  const auto potential_result_edge = mergeable_edges[0];
  const auto [source_id, sink_id] = mSourceAndSinkIds;
  if (potential_result_edge.DstId() == source_id || potential_result_edge.DstId() == sink_id) {
    return std::nullopt;
  }

  // Check if src --> abc_bdy is a potential buddy edge
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (!IsPotentialBuddyEdge(src, potential_result_edge)) return std::nullopt;

  const auto opp_dir_edges = src.FindEdgesInDirection(Kmer::RevOrdering(ord));
  // NOLINTBEGIN(readability-braces-around-statements)
  if (opp_dir_edges.empty()) return potential_result_edge;
  if (opp_dir_edges.size() > 1) return std::nullopt;
  // NOLINTEND(readability-braces-around-statements)

  // Check if src --> abc_bdy is a potential buddy edge
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (!IsPotentialBuddyEdge(src, opp_dir_edges[0])) return std::nullopt;

  return potential_result_edge;
}

auto Graph::IsPotentialBuddyEdge(const Node& src, const Edge& conn) const -> bool {
  // conn is an outgoing edge from src node. conn's dst node is called nbour for "neighbour" node.
  // * nbour node has atmost 2 and atleast 1 outgoing edges and not contain any self loops.
  // * One of nbour outgoing edges must be a mirror of src --> nbour. i.e. nbour --> src.
  // * Another nbour outgoing edge if present, must in the opposite direction of nbour --> src to nbour's nbour `nnb`
  // * Neighbour's neighbour `nnb` if present, can atmost have 2 outgoing edges
  // If all of these checks pass, then `conn` is a potential buddy edge from src --> neighbour
  const auto nbour_itr = mNodes.find(conn.DstId());
  LANCET_ASSERT(nbour_itr != mNodes.end())
  LANCET_ASSERT(nbour_itr->second != nullptr)
  const Node& nbour = *nbour_itr->second;

  // Check edge case where the only nodes between src and nbour are each other
  if (src.NumOutEdges() == 1 && nbour.NumOutEdges() == 1) {
    const auto edge_from_src = std::vector<Edge>(src.cbegin(), src.cend())[0];
    const auto edge_from_nbour = std::vector<Edge>(nbour.cbegin(), nbour.cend())[0];
    if (edge_from_src.DstId() == nbour.Identifier() && edge_from_nbour.DstId() == src.Identifier()) {
      return false;
    }
  }

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (nbour.NumOutEdges() > 2 || nbour.NumOutEdges() == 0 || nbour.HasSelfLoop()) return false;

  const auto expected_nbour2src = conn.MirrorEdge();
  const auto start_sign_nbour2src = expected_nbour2src.SrcSign();
  const auto dir_nbour2src = start_sign_nbour2src == nbour.SignFor(Kmer::Ordering::DEFAULT) ? Kmer::Ordering::DEFAULT
                                                                                            : Kmer::Ordering::OPPOSITE;
  const auto nb_edges_in_nbour2src_dir = nbour.FindEdgesInDirection(dir_nbour2src);
  if (nb_edges_in_nbour2src_dir.size() != 1 || nb_edges_in_nbour2src_dir[0] != expected_nbour2src) {
    return false;
  }

  const auto nb_edges_in_opp_dir = nbour.FindEdgesInDirection(Kmer::RevOrdering(dir_nbour2src));
  // Check if nbour loops back in a cycle to src node in opposite direction again
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (nb_edges_in_opp_dir.size() != 1 || nb_edges_in_opp_dir[0].DstId() == conn.SrcId()) return false;

  const auto nnb_itr = mNodes.find(nb_edges_in_opp_dir[0].DstId());
  LANCET_ASSERT(nnb_itr != mNodes.end())
  LANCET_ASSERT(nnb_itr->second != nullptr)
  return nnb_itr->second->NumOutEdges() <= 2;
}

void Graph::RemoveTips(const usize component_id) {
  usize total_tips = 0;
  usize curr_tips = 1;

  std::vector<NodeID> remove_nids;
  remove_nids.reserve(mNodes.size());
  // remove tips and compress at least once. compression after tip removal
  // can produce new tips in the graph, so recursively remove tips from
  // the graph until there are no longer any tips left
  while (curr_tips > 0) {
    remove_nids.clear();

    std::ranges::for_each(mNodes, [&remove_nids, &component_id, this](NodeTable::const_reference item) {
      const auto [source_id, sink_id] = this->mSourceAndSinkIds;
      // NOLINTBEGIN(readability-braces-around-statements)
      if (item.second->GetComponentId() != component_id || item.second->NumOutEdges() > 1) return;
      if (item.first == source_id || item.first == sink_id) return;
      // NOLINTEND(readability-braces-around-statements)

      const auto uniq_seq_len = item.second->SeqLength() - mCurrK + 1;
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (uniq_seq_len >= this->mCurrK) return;

      remove_nids.emplace_back(item.first);
    });

    if (!remove_nids.empty()) {
      total_tips += curr_tips;
      RemoveNodes(absl::MakeConstSpan(remove_nids));
      CompressGraph(component_id);
    }

    curr_tips = remove_nids.size();
  }

  if (total_tips > 0) {
    // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
    const auto region_str = mRegion->ToSamtoolsRegion();
    LOG_TRACE("Removed {} tips for {} in comp{} with k={}", total_tips, region_str, component_id, mCurrK)
  }
}

auto Graph::FindSource(const usize component_id) const -> RefAnchor {
  RefAnchor result{.mAnchorId = 0, .mRefOffset = 0, .mFoundAnchor = false};

  for (usize ref_idx = 0; ref_idx < mRefNodeIds.size(); ++ref_idx) {
    const auto itr = mNodes.find(mRefNodeIds[ref_idx]);
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (itr == mNodes.end()) continue;

    LANCET_ASSERT(itr->second != nullptr)
    if (itr->second->GetComponentId() != component_id || itr->second->TotalReadSupport() < mParams.mMinAnchorCov) {
      continue;
    }

    result.mAnchorId = itr->first;
    result.mRefOffset = ref_idx;
    result.mFoundAnchor = true;
    break;
  }

  return result;
}

auto Graph::FindSink(const usize component_id) const -> RefAnchor {
  RefAnchor result{.mAnchorId = 0, .mRefOffset = 0, .mFoundAnchor = false};

  for (i64 ref_idx = static_cast<i64>(mRefNodeIds.size() - 1); ref_idx >= 0; --ref_idx) {
    const auto itr = mNodes.find(mRefNodeIds[ref_idx]);
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (itr == mNodes.end()) continue;

    LANCET_ASSERT(itr->second != nullptr)
    if (itr->second->GetComponentId() != component_id || itr->second->TotalReadSupport() < mParams.mMinAnchorCov) {
      continue;
    }

    result.mAnchorId = itr->first;
    result.mRefOffset = static_cast<usize>(ref_idx);
    result.mFoundAnchor = true;
    break;
  }

  return result;
}

// ============================================================================
//  HasCycle — O(V+E) Three-Color DFS on the Flat Adjacency List
// ============================================================================
//
// ALGORITHM
// ---------
// Standard directed-graph cycle detection using three colors:
//   WHITE (0) = unvisited
//    GRAY (1) = on the current DFS stack (ancestor in the current path)
//   BLACK (2) = fully explored (all descendants visited)
//
// A cycle exists iff DFS encounters a "back edge" — an edge leading to a
// GRAY state. GRAY states are ancestors on the current DFS path, so a back
// edge forms a cycle.
//
// WHY THIS REPLACES THE OLD APPROACH
// ------------------------------------
// The old HasCycle used backtracking (erase from visited set on return),
// which explored ALL paths from source — exponential in high-branching
// graphs. Three-color DFS visits each state exactly once: O(V+E).
// Profile data showed ~51.6s in the old HasCycle; this should be <1ms.
//
// BIDIRECTED SIGN CONTINUITY
// ----------------------------
// In the BCALM2 bidirected model, walks must satisfy sign continuity:
//   edge.DstSign must match the next edge.SrcSign
//
// We track state as (node_flat_idx, sign), so a node visited via '+' and
// via '-' are different DFS states. The TraversalIndex adjacency list is
// already partitioned by (node, sign), so edge iteration naturally respects
// sign continuity.
//
//   Example: DFS from source(+)
//   ┌──────────────────────────────────────────────────────┐
//   │  Visit state (A,+) → GRAY                           │
//   │    Edge to (B,+) → WHITE → push (B,+)               │
//   │      Edge to (C,-) → WHITE → push (C,-)             │
//   │        Edge to (A,+) → GRAY! → CYCLE FOUND          │
//   │      Edge to (A,-) → WHITE → push (A,-)             │ ← same node, different sign
//   │        ...no back edge → BLACK                       │
//   └──────────────────────────────────────────────────────┘
//
// FUTURE: Tarjan's SCC could provide cycle sizes and weakest edges for
// smarter retry-vs-skip decisions. For now, just detect presence/absence.
//
auto Graph::HasCycle(const TraversalIndex& idx) const -> bool {
  enum Color : u8 { WHITE = 0, GRAY = 1, BLACK = 2 };

  // Flat color array indexed by state_idx = node_flat_idx * 2 + sign_offset.
  // O(1) access per state — no hash lookups, cache-line-friendly.
  std::vector<u8> color(idx.NumStates(), WHITE);

  // Iterative DFS stack frame: current state + position in its outgoing edge list.
  // Using an explicit stack avoids recursion depth limits on large graphs.
  struct DfsFrame {
    u32 mStateIdx;
    u32 mEdgePos;  // next edge to explore within this state's adjacency range
  };

  std::vector<DfsFrame> stack;
  stack.reserve(idx.NumNodes());

  color[idx.mSrcState] = GRAY;
  stack.push_back({idx.mSrcState, 0});

  while (!stack.empty()) {
    auto& frame = stack.back();
    const auto& range = idx.mAdjRanges[frame.mStateIdx];

    // All children of this state explored → mark BLACK and backtrack
    if (frame.mEdgePos >= range.mCount) {
      color[frame.mStateIdx] = BLACK;
      stack.pop_back();
      continue;
    }

    // Examine the next outgoing edge from this state
    const auto& out = idx.mAdjList[range.mStart + frame.mEdgePos];
    frame.mEdgePos++;

    if (color[out.mDstState] == GRAY) return true;  // Back edge → cycle!
    if (color[out.mDstState] != WHITE) continue;     // BLACK → already finished, skip

    // WHITE → unvisited. Push child onto DFS stack.
    color[out.mDstState] = GRAY;
    stack.push_back({out.mDstState, 0});
  }

  return false;
}

// ============================================================================
//  BuildTraversalIndex — Construct Flat CSR Adjacency List
// ============================================================================
//
// Converts the hash-map-based NodeTable into a contiguous, integer-indexed
// adjacency list for a single connected component. This is built ONCE after
// all graph mutations (pruning, compression, tip removal) are complete.
//
// CONSTRUCTION PHASES
// --------------------
//  Phase 1: Assign contiguous u32 indices to nodes in this component.
//  Phase 2: Count outgoing edges per state (for CSR range sizing).
//  Phase 3: Compute prefix-sum offsets for each state's edge range.
//  Phase 4: Fill the adjacency list and assign edge ordinals.
//  Phase 5: Set source and sink state indices.
//
// COST: O(V + E) time and memory. The nid_to_flat hash map is only used
// during construction; all subsequent operations are flat-array-only.
//
auto Graph::BuildTraversalIndex(const usize component_id) const -> TraversalIndex {
  TraversalIndex idx;

  // Phase 1: Assign contiguous u32 indices to nodes in this component
  absl::flat_hash_map<NodeID, u32> nid_to_flat;
  nid_to_flat.reserve(mNodes.size());

  for (const auto& [nid, node_ptr] : mNodes) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (node_ptr->GetComponentId() != component_id) continue;
    const auto flat = static_cast<u32>(idx.mNodes.size());
    idx.mNodes.push_back(node_ptr.get());
    idx.mNodeIds.push_back(nid);
    nid_to_flat.emplace(nid, flat);
  }

  const u32 num_nodes = idx.NumNodes();
  const u32 num_states = num_nodes * 2;
  idx.mAdjRanges.resize(num_states, {0, 0});

  // Phase 2: Count outgoing edges per state (for range sizing)
  for (u32 ni = 0; ni < num_nodes; ni++) {
    const Node* node = idx.mNodes[ni];
    for (const Edge& edge : *node) {
      // Only count edges whose destination is in this component
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (nid_to_flat.find(edge.DstId()) == nid_to_flat.end()) continue;
      const u32 state = TraversalIndex::MakeState(ni, edge.SrcSign());
      idx.mAdjRanges[state].mCount++;
    }
  }

  // Phase 3: Compute starting offsets via prefix sum
  u32 offset = 0;
  for (u32 s = 0; s < num_states; s++) {
    idx.mAdjRanges[s].mStart = offset;
    offset += idx.mAdjRanges[s].mCount;
    idx.mAdjRanges[s].mCount = 0;  // reset count; re-filled in phase 4
  }
  idx.mAdjList.resize(offset);

  // Phase 4: Fill adjacency list entries and assign edge ordinals
  absl::flat_hash_map<Edge, u32> edge_to_ordinal;
  edge_to_ordinal.reserve(offset);

  for (u32 ni = 0; ni < num_nodes; ni++) {
    const Node* node = idx.mNodes[ni];
    for (const Edge& edge : *node) {
      const auto dst_it = nid_to_flat.find(edge.DstId());
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (dst_it == nid_to_flat.end()) continue;

      const u32 src_state = TraversalIndex::MakeState(ni, edge.SrcSign());
      const u32 dst_state = TraversalIndex::MakeState(dst_it->second, edge.DstSign());

      // Assign or reuse edge ordinal (edges appear at both endpoints as forward + mirror)
      u32 ordinal = 0;
      auto [it, inserted] = edge_to_ordinal.emplace(edge, static_cast<u32>(idx.mOrigEdges.size()));
      if (inserted) {
        ordinal = it->second;
        idx.mOrigEdges.push_back(edge);
      } else {
        ordinal = it->second;
      }

      auto& range = idx.mAdjRanges[src_state];
      idx.mAdjList[range.mStart + range.mCount] = {dst_state, ordinal};
      range.mCount++;
    }
  }

  // Phase 5: Set source and sink states
  const auto [source_id, sink_id] = mSourceAndSinkIds;
  const auto src_flat = nid_to_flat.at(source_id);
  const auto snk_flat = nid_to_flat.at(sink_id);
  const auto src_sign = idx.mNodes[src_flat]->SignFor(Kmer::Ordering::DEFAULT);
  idx.mSrcState = TraversalIndex::MakeState(src_flat, src_sign);
  idx.mSnkNodeIdx = snk_flat;

  return idx;
}

// ============================================================================
//  ComputeGraphComplexity — O(V+E) Metrics for Debug Logging
// ============================================================================
//
// Computes lightweight graph topology metrics that correlate with runtime:
//
//  Cyclomatic complexity (M = E - V + 1): number of independent cycles.
//    M=0 → linear chain, M=1 → single variant bubble, M>>1 → STR hairball.
//
//  Edge-to-node density (E/V): a clean graph has E/V ≈ 1.0. Values >1.5
//    indicate branching that causes path enumeration explosion.
//
//  Max single-direction degree: maximum outgoing edges in any one sign
//    direction. Hub nodes with high degree are direct BFS blowup predictors.
//
//  Branch points: nodes with ≥2 outgoing edges in some direction. These
//    are the "decision points" that cause combinatorial path blowup.
//
auto Graph::ComputeGraphComplexity(const usize component_id) const -> GraphComplexity {
  GraphComplexity cx;
  for (const auto& [nid, node_ptr] : mNodes) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (node_ptr->GetComponentId() != component_id) continue;
    cx.mNumNodes++;

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
    cx.mNumEdges += dflt_dir_edges + oppo_dir_edges;
    const auto max_dir = std::max(dflt_dir_edges, oppo_dir_edges);
    cx.mMaxSingleDirDegree = std::max(cx.mMaxSingleDirDegree, max_dir);
    if (dflt_dir_edges >= 2 || oppo_dir_edges >= 2) cx.mNumBranchPoints++;
  }

  // Each edge stored at both endpoints (forward + mirror) → divide by 2
  cx.mNumEdges /= 2;
  cx.mCyclomaticComplexity =
      (cx.mNumEdges >= cx.mNumNodes) ? (cx.mNumEdges - cx.mNumNodes + 1) : 0;
  cx.mEdgeToNodeDensity = cx.mNumNodes > 0
      ? static_cast<f64>(cx.mNumEdges) / static_cast<f64>(cx.mNumNodes) : 0.0;
  return cx;
}


auto Graph::MarkConnectedComponents() -> std::vector<ComponentInfo> {
  usize current_component = 0;
  std::vector<ComponentInfo> results_info;

#ifdef LANCET_DEVELOP_MODE
  static const auto is_unassigned = [](NodeTable::const_reference item) { return item.second->GetComponentId() == 0; };
#endif

  // Check that all nodes are component zero before we start
  LANCET_ASSERT(static_cast<usize>(std::ranges::count_if(mNodes, is_unassigned)) == mNodes.size())

  for (NodeTable::reference item : mNodes) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (item.second->GetComponentId() != 0) continue;

    current_component++;
    results_info.emplace_back(ComponentInfo{.mCompId = current_component, .mNumNodes = 0});

    absl::chunked_queue<Node*, 128, 1024> connected_nodes;
    connected_nodes.push_back(item.second.get());

    while (!connected_nodes.empty()) {
      auto* current_node = connected_nodes.front();
      LANCET_ASSERT(current_node != nullptr)

      if (current_node->GetComponentId() != 0) {
        connected_nodes.pop_front();
        continue;
      }

      current_node->SetComponentId(current_component);
      results_info[current_component - 1].mNumNodes += 1;
      for (const Edge& edge : *current_node) {
        const auto neighbour_itr = mNodes.find(edge.DstId());
        LANCET_ASSERT(neighbour_itr != mNodes.end())
        LANCET_ASSERT(neighbour_itr->second != nullptr)
        connected_nodes.push_back(neighbour_itr->second.get());
      }

      connected_nodes.pop_front();
    }
  }

  const auto total_num_nodes = static_cast<f64>(mNodes.size());
  std::ranges::for_each(results_info, [&total_num_nodes](ComponentInfo& cinfo) {
    cinfo.mPctNodes = 100.0 * (static_cast<f64>(cinfo.mNumNodes) / total_num_nodes);
  });

  std::ranges::sort(results_info, [](const ComponentInfo& lhs, const ComponentInfo& rhs) -> bool {
    return lhs.mNumNodes > rhs.mNumNodes;
  });

  // Check that none of the nodes are component zero after we are done
  LANCET_ASSERT(static_cast<usize>(std::ranges::count_if(mNodes, is_unassigned)) == 0)
  return results_info;
}

void Graph::RemoveLowCovNodes(const usize component_id) {
  std::vector<NodeID> remove_nids;
  remove_nids.reserve(mNodes.size());

  std::ranges::for_each(std::as_const(mNodes), [&remove_nids, &component_id, this](NodeTable::const_reference item) {
    const auto [source_id, sink_id] = this->mSourceAndSinkIds;
    // NOLINTBEGIN(readability-braces-around-statements)
    if (item.second->GetComponentId() != component_id) return;
    if (item.first == source_id || item.first == sink_id) return;
    // NOLINTEND(readability-braces-around-statements)

    const auto is_nml_singleton = item.second->NormalReadSupport() == 1;
    const auto is_tmr_singleton = item.second->TumorReadSupport() == 1;
    const auto total_sample_cov = item.second->TotalReadSupport();

    if ((is_nml_singleton && is_tmr_singleton) || total_sample_cov < this->mParams.mMinNodeCov) {
      remove_nids.emplace_back(item.first);
    }
  });

  if (!remove_nids.empty()) {
    // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
    const auto region_str = mRegion->ToSamtoolsRegion();
    LOG_TRACE("Removing {:.4f}% (or) {} low cov nodes for {} in comp{} with k={}",
              100.0 * (static_cast<f64>(remove_nids.size()) / static_cast<f64>(mNodes.size())), remove_nids.size(),
              region_str, component_id, mCurrK)

    RemoveNodes(absl::MakeConstSpan(remove_nids));
  }
}

void Graph::RemoveNode(NodeTable::iterator itr) {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (itr == mNodes.end()) return;

  // remove all incoming edges to the node first
  for (const Edge& conn : *itr->second) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (conn.IsSelfLoop()) continue;

    auto nbour_itr = mNodes.find(conn.DstId());
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (nbour_itr != mNodes.end()) nbour_itr->second->EraseEdge(conn.MirrorEdge());
  }

  mNodes.erase(itr);
}

void Graph::RemoveNodes(absl::Span<const NodeID> node_ids) {
  std::ranges::for_each(node_ids, [this](const NodeID nid) { this->RemoveNode(this->mNodes.find(nid)); });
}

void Graph::BuildGraph(absl::flat_hash_set<MateMer>& mate_mers) {
  mRefNodeIds.clear();
  const auto ref_nodes = AddNodes(mRegion->SeqView(), Label::REFERENCE);
  mRefNodeIds.reserve(ref_nodes.size());
  std::ranges::transform(ref_nodes, std::back_inserter(mRefNodeIds),
                         [](const Node* node) -> NodeID { return node->Identifier(); });

  // Expected errors = floor(err_sum) https://www.drive5.com/usearch/manual/exp_errs.html
  // See https://doi.org/10.1093/bioinformatics/btv401 for proof on expected errors
  // Add support for only high quality kmers. If the expected error is > MAX_AGG_ERR,
  // then skip adding any read support for those kmers, so they can be removed later
  static const auto error_summer = [](const f64 sum, const u8 bql) -> f64 { return sum + hts::PhredToErrorProb(bql); };
  static const auto is_low_qual_kmer = [](absl::Span<const u8> quals) -> bool {
    const f64 expected_errors = std::floor(std::accumulate(quals.cbegin(), quals.cend(), 0.0, error_summer));
    return static_cast<i64>(expected_errors) > 0;
  };

  mate_mers.clear();
  for (const auto& read : mReads) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (!read.PassesAlnFilters()) continue;

    usize offset = 0;
    auto added_nodes = AddNodes(read.SeqView(), read.SrcLabel());

    std::ranges::for_each(added_nodes, [&read, &offset, &mate_mers, this](Node* node) {
      MateMer mm_pair{.mQname = read.QnameView(), .mTagKind = read.TagKind(), .mKmerHash = node->Identifier()};
      const auto curr_qual = read.QualView().subspan(offset, this->mCurrK);
      offset++;

      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (is_low_qual_kmer(curr_qual) || mate_mers.contains(mm_pair)) return;
      node->IncrementReadSupport(read.SrcLabel());
      mate_mers.emplace(std::move(mm_pair));
    });
  }
}

auto Graph::AddNodes(std::string_view sequence, const Label label) -> std::vector<Node*> {
  std::vector<Node*> result;
  const auto kplus_ones = SlidingView(sequence, mCurrK + 1);
  result.reserve(kplus_ones.size() + 1);

  for (usize mer_idx = 0; mer_idx < kplus_ones.size(); ++mer_idx) {
    const auto seq1 = absl::ClippedSubstr(kplus_ones[mer_idx], 0, mCurrK);
    const auto seq2 = absl::ClippedSubstr(kplus_ones[mer_idx], 1, mCurrK);

    auto left_mer = Kmer(seq1);
    auto right_mer = Kmer(seq2);
    const auto left_id = left_mer.Identifier();
    const auto right_id = right_mer.Identifier();

    mNodes.try_emplace(left_id, std::make_unique<Node>(std::move(left_mer), label));
    mNodes.try_emplace(right_id, std::make_unique<Node>(std::move(right_mer), label));

    auto& first = mNodes.at(left_id);
    auto& second = mNodes.at(right_id);

    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (mer_idx == 0) result.emplace_back(first.get());

    static constexpr auto dflt_order = Kmer::Ordering::DEFAULT;
    const auto fwd_edge = MakeFwdEdgeKind({first->SignFor(dflt_order), second->SignFor(dflt_order)});
    first->EmplaceEdge(NodeIDPair{left_id, right_id}, fwd_edge);
    second->EmplaceEdge(NodeIDPair{right_id, left_id}, RevEdgeKind(fwd_edge));

    result.emplace_back(second.get());
  }

  return result;
}

auto Graph::HasExactOrApproxRepeat(std::string_view seq, usize window) -> bool {
  const auto klen_seqs = SlidingView(seq, window);
  static constexpr usize NUM_ALLOWED_MISMATCHES = 3;

  return HasExactRepeat(absl::MakeConstSpan(klen_seqs)) ||
         HasApproximateRepeat(absl::MakeConstSpan(klen_seqs), NUM_ALLOWED_MISMATCHES);
}

auto Graph::RefAnchorLength(const RefAnchor& source, const RefAnchor& sink, usize currk) -> usize {
  return sink.mRefOffset - source.mRefOffset + currk;
}

#ifdef LANCET_DEVELOP_MODE
auto Graph::ToString(const State state) -> std::string {
  switch (state) {
    case FIRST_LOW_COV_REMOVAL:
      return "low_cov_removal1";
    case FOUND_REF_ANCHORS:
      return "found_ref_anchors";
    case FIRST_COMPRESSION:
      return "compression1";
    case SECOND_LOW_COV_REMOVAL:
      return "low_cov_removal2";
    case SECOND_COMPRESSION:
      return "compression2";
    case SHORT_TIP_REMOVAL:
      return "short_tip_removal";
    default:
      break;
  }

  return "fully_pruned";
}
#endif

void Graph::WriteDot([[maybe_unused]] State state, usize comp_id) {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mParams.mOutGraphsDir.empty()) return;

#ifdef LANCET_DEVELOP_MODE
  const auto graph_state = ToString(state);
#else
  const auto* graph_state = "fully_pruned";
#endif

  using namespace std::string_view_literals;
  const auto win_id = fmt::format("{}_{}_{}", mRegion->ChromName(), mRegion->StartPos1(), mRegion->EndPos1());
  const auto fname = fmt::format("dbg__{}__{}__k{}__comp{}.dot", win_id, graph_state, mCurrK, comp_id);

  const auto out_path = mParams.mOutGraphsDir / "dbg_graph" / fname;
  std::filesystem::create_directories(mParams.mOutGraphsDir / "dbg_graph");
  SerializeToDot(mNodes, out_path, comp_id, {mSourceAndSinkIds.cbegin(), mSourceAndSinkIds.cend()});
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Graph::SerializeToDot(const NodeTable& graph, const std::filesystem::path& out_path, const usize comp_id,
                           const NodeIdSet& nodes_highlight, const EdgeSet& edges_highlight,
                           const NodeIdSet& nodes_background, const EdgeSet& edges_background) {
  std::ofstream out_handle(out_path, std::ios::trunc);
  using namespace std::string_view_literals;

  out_handle << R"raw(strict digraph G {
graph [layout=neato,bgcolor=black,size="120,180",ratio=compress,rankdir=LR,overlap=vpsc,overlap_shrink=true,start=self];
node [style=filled,fontsize=2,width=2,height=2,fixedsize=false];
edge [color=gray,fontsize=8,fontcolor=floralwhite,len=3,fixedsize=false,headclip=true,tailclip=true];
)raw"sv;

  fmt::print(out_handle, "subgraph {} {{\n", out_path.stem().string());

  for (NodeTable::const_reference item : graph) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (item.second->GetComponentId() != comp_id) continue;

    const auto dflt_seq = item.second->SequenceFor(Kmer::Ordering::DEFAULT);
    const auto oppo_seq = item.second->SequenceFor(Kmer::Ordering::OPPOSITE);
    const auto rev_oppo_seq = std::string(oppo_seq.crbegin(), oppo_seq.crend());
    const auto sign_dflt = item.second->SignFor(Kmer::Ordering::DEFAULT) == Kmer::Sign::PLUS ? '+' : '-';
    const auto is_background_node = nodes_background.contains(item.first);
    // NOLINTBEGIN(readability-avoid-nested-conditional-operator)
    const auto fill_color = is_background_node                     ? "darkgray"sv
                            : nodes_highlight.contains(item.first) ? "orchid"sv
                            : item.second->IsShared()              ? "steelblue"sv
                            : item.second->IsTumorOnly()           ? "indianred"sv
                            : item.second->IsNormalOnly()          ? "mediumseagreen"sv
                                                                   : "lightblue"sv;
    // NOLINTEND(readability-avoid-nested-conditional-operator)

    fmt::print(out_handle, R"raw({} [shape=circle fillcolor={} label="{}\n{}\n {}:{}\nlength={}\ncoverage={}"]
)raw",
               item.first, fill_color, dflt_seq, rev_oppo_seq, item.first, sign_dflt, item.second->Length(),
               item.second->TotalReadSupport());

    for (const Edge& conn : *item.second) {
      const auto src_sign = conn.SrcSign() == Kmer::Sign::PLUS ? '+' : '-';
      const auto dst_sign = conn.DstSign() == Kmer::Sign::PLUS ? '+' : '-';
      const auto is_background_edge = edges_background.contains(conn);
      const auto is_highlight_edge = edges_highlight.contains(conn);
      fmt::print(out_handle, R"raw({} -> {} [taillabel="{}" headlabel="{}" style="{}"{}]
)raw",
                 conn.SrcId(), conn.DstId(), src_sign, dst_sign, is_background_edge ? "dotted"sv : "solid"sv,
                 is_highlight_edge ? R"raw( color="goldenrod")raw"sv : ""sv);
    }
  }

  out_handle << "}\n}\n"sv;
  out_handle.close();
}

}  // namespace lancet::cbdg
