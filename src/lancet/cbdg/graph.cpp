#include "lancet/cbdg/graph.h"

#include <algorithm>
#include <deque>
#include <fstream>

#include "absl/hash/hash.h"
#include "absl/strings/str_join.h"
#include "absl/strings/string_view.h"
#include "lancet/base/assert.h"
#include "lancet/base/hash.h"
#include "lancet/base/logging.h"
#include "lancet/base/repeat.h"
#include "lancet/base/rev_comp.h"
#include "lancet/base/sliding.h"
#include "lancet/base/timer.h"
#include "lancet/cbdg/max_flow.h"
#include "spdlog/fmt/fmt.h"
#include "spdlog/fmt/ostr.h"

namespace lancet::cbdg {

/// https://github.com/GATB/bcalm/blob/v2.2.3/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
auto Graph::BuildComponentHaplotypes(RegionPtr region, ReadList reads) -> Result {
  mReads = reads;
  mRegion = std::move(region);

  Timer timer;
  GraphHaps per_comp_haplotypes;
  std::string_view ref_anchor_seq;
  std::vector<usize> anchor_start_idxs;
  absl::flat_hash_set<MateMer> mate_mers;

  static constexpr usize DEFAULT_MIN_ANCHOR_LENGTH = 150;
  static constexpr f64 DEFAULT_PCT_NODES_NEEDED = 10.0;

  const auto reg_str = mRegion->ToSamtoolsRegion();
  mCurrK = mParams.mMinKmerLen - 2;

IncrementKmerAndRetry:
  while (per_comp_haplotypes.empty() && mCurrK < mParams.mMaxKmerLen) {
    mCurrK += 2;
    timer.Reset();
    mAverageCov = 0.0;
    mSourceAndSinkIds = {0, 0};

    // NOLINTNEXTLINE(readability-braces-around-statements,cppcoreguidelines-avoid-goto)
    if (HasExactOrApproxRepeat(mRegion->SeqView(), mCurrK)) goto IncrementKmerAndRetry;

    BuildGraph(mate_mers);
    LOG_TRACE("Built graph for {} with k={}, nodes={}, reads={}", reg_str, mCurrK, mNodes.size(), mReads.size())

    RemoveLowCovNodes(0);
    mNodes.rehash(0);
    WriteDotDevelop(FIRST_LOW_COV_REMOVAL, 0);

    const auto components = MarkConnectedComponents();
    per_comp_haplotypes.reserve(components.size());
    anchor_start_idxs.reserve(components.size());
    LOG_TRACE("Found {} connected components in graph for {} with k={}", components.size(), reg_str, mCurrK)

    for (const auto& cinfo : components) {
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (cinfo.mPctNodes < DEFAULT_PCT_NODES_NEEDED) continue;

      const auto comp_id = cinfo.mCompId;
      const auto source = FindSource(comp_id);
      const auto sink = FindSink(comp_id);

      if (!source.mFoundAnchor || !sink.mFoundAnchor || source.mAnchorId == sink.mAnchorId) {
        LOG_TRACE("Skipping comp{} in graph for {} because source and sink were not found", comp_id, reg_str)
        continue;
      }

      const auto current_anchor_length = RefAnchorLength(source, sink, mCurrK);
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (current_anchor_length < DEFAULT_MIN_ANCHOR_LENGTH) continue;

      LOG_TRACE("Found {}bp anchor for {} comp={} with k={}", current_anchor_length, reg_str, comp_id, mCurrK)

      std::vector<std::string> haplotypes;
      mSourceAndSinkIds = NodeIDPair{source.mAnchorId, sink.mAnchorId};
      ref_anchor_seq = mRegion->SeqView().substr(source.mRefOffset, current_anchor_length);
      WriteDotDevelop(FOUND_REF_ANCHORS, comp_id);

      if (HasCycle()) {
        LOG_TRACE("Graph cycle found for {} comp={} with k={}", reg_str, comp_id, mCurrK)
        goto IncrementKmerAndRetry;  // NOLINT(cppcoreguidelines-avoid-goto)
      }

      CompressGraph(comp_id);
      WriteDotDevelop(FIRST_COMPRESSION, comp_id);
      RemoveLowCovNodes(comp_id);
      WriteDotDevelop(SECOND_LOW_COV_REMOVAL, comp_id);
      CompressGraph(comp_id);
      WriteDotDevelop(SECOND_COMPRESSION, comp_id);
      RemoveTips(comp_id);
      WriteDotDevelop(SHORT_TIP_REMOVAL, comp_id);

      if (HasCycle()) {
        LOG_TRACE("Graph cycle found for {} comp={} with k={}", reg_str, comp_id, mCurrK)
        goto IncrementKmerAndRetry;  // NOLINT(cppcoreguidelines-avoid-goto)
      }

      WriteDot(State::FULLY_PRUNED_GRAPH, comp_id);
      LOG_TRACE("Starting Edmond Karp traversal for {} with k={}, num_nodes={}", reg_str, mCurrK, mNodes.size())
      MaxFlow max_flow(&mNodes, mSourceAndSinkIds, mCurrK);
      auto path_seq = max_flow.NextPath();

      while (path_seq) {
        LOG_TRACE("Assembled {}bp path sequence for {} with k={}", path_seq->length(), reg_str, mCurrK)
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
  }

  static const auto summer = [](const u64 sum, const auto& comp_haps) -> u64 { return sum + comp_haps.size() - 1; };
  const auto num_asm_haps = std::accumulate(per_comp_haplotypes.cbegin(), per_comp_haplotypes.cend(), 0, summer);
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

  std::vector<NodeID> nids_to_remove;
  nids_to_remove.reserve(mNodes.size());
  // remove tips and compress at least once. compression after tip removal
  // can produce new tips in the graph, so recursively remove tips from
  // the graph until there are no longer any tips left
  while (curr_tips > 0) {
    nids_to_remove.clear();

    std::ranges::for_each(mNodes, [&nids_to_remove, &component_id, this](NodeTable::const_reference item) {
      const auto [source_id, sink_id] = this->mSourceAndSinkIds;
      // NOLINTBEGIN(readability-braces-around-statements)
      if (item.second->GetComponentId() != component_id || item.second->NumOutEdges() > 1) return;
      if (item.first == source_id || item.first == sink_id) return;
      // NOLINTEND(readability-braces-around-statements)

      const auto uniq_seq_len = item.second->SeqLength() - mCurrK + 1;
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (uniq_seq_len >= this->mCurrK) return;

      nids_to_remove.emplace_back(item.first);
    });

    if (!nids_to_remove.empty()) {
      total_tips += curr_tips;
      RemoveNodes(absl::MakeConstSpan(nids_to_remove));
      CompressGraph(component_id);
    }

    curr_tips = nids_to_remove.size();
  }

  if (total_tips > 0) {
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

auto Graph::HasCycle() const -> bool {
  const auto src_itr = mNodes.find(mSourceAndSinkIds[0]);
  LANCET_ASSERT(src_itr != mNodes.end())
  LANCET_ASSERT(src_itr->second != nullptr)

  bool cycle_found = false;
  usize recursion_count = 0;
  absl::flat_hash_set<NodeID> traversed;
  traversed.reserve(mNodes.size());

  HasCycle(*src_itr->second, traversed, cycle_found, recursion_count);
  return cycle_found;
}

void Graph::HasCycle(const Node& node, NodeIdSet& traversed, bool& found_cycle, usize& recursion_depth) const {
  const auto node_default_sign = node.SignFor(Kmer::Ordering::DEFAULT);
  traversed.insert(node.Identifier());

  // If we get here, then the graph is too complex and we are
  // stuck in a deep recursion context with no way out.
  const auto max_recursion_limit = mNodes.size() * mNodes.size();
  if (recursion_depth > max_recursion_limit) {
    found_cycle = true;
    return;
  }

  for (const Edge& conn : node) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (conn.SrcSign() != node_default_sign) continue;
    if (traversed.find(conn.DstId()) != traversed.end()) {
      found_cycle = true;
      return;
    }

    const auto neighbour_itr = mNodes.find(conn.DstId());
    LANCET_ASSERT(neighbour_itr != mNodes.end())
    LANCET_ASSERT(neighbour_itr->second != nullptr)

    recursion_depth++;
    HasCycle(*neighbour_itr->second, traversed, found_cycle, recursion_depth);
  }

  traversed.erase(node.Identifier());
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

    std::deque<Node*> connected_nodes;
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
  // min_node_cov -> minimum coverage required for each node.
  // min_ratio_cov -> combined sample coverage * MIN_NODE_COV_RATIO for each node
  const auto min_ratio_cov = static_cast<u32>(std::floor(mParams.mMinNodeCovRatio * mAverageCov));
  const auto min_req_cov = std::max(mParams.mMinNodeCov, min_ratio_cov);

  std::vector<NodeID> nodes_to_remove;
  nodes_to_remove.reserve(mNodes.size());

  std::ranges::for_each(std::as_const(mNodes),
                        [&nodes_to_remove, &component_id, &min_req_cov, this](NodeTable::const_reference item) {
                          const auto [source_id, sink_id] = this->mSourceAndSinkIds;
                          // NOLINTBEGIN(readability-braces-around-statements)
                          if (item.second->GetComponentId() != component_id) return;
                          if (item.first == source_id || item.first == sink_id) return;
                          // NOLINTEND(readability-braces-around-statements)

                          const auto is_nml_singleton = item.second->NormalReadSupport() == 1;
                          const auto is_tmr_singleton = item.second->TumorReadSupport() == 1;
                          const auto total_sample_cov = item.second->TotalReadSupport();

                          if ((is_nml_singleton && is_tmr_singleton) || total_sample_cov < min_req_cov) {
                            nodes_to_remove.emplace_back(item.first);
                          }
                        });

  if (!nodes_to_remove.empty()) {
    const auto region_str = mRegion->ToSamtoolsRegion();
    LOG_TRACE("Removing {:.4f}% (or) {} low coverage nodes for {} in comp{} with k={}",
              100.0 * (static_cast<f64>(nodes_to_remove.size()) / static_cast<f64>(mNodes.size())),
              nodes_to_remove.size(), region_str, component_id, mCurrK)

    RemoveNodes(absl::MakeConstSpan(nodes_to_remove));
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
  usize nsample_bases = 0;
  usize max_num_kmers = 0;
  std::vector<Label> labels;
  std::vector<SeqMers> kplus_ones;
  labels.reserve(mReads.size() + 1);
  kplus_ones.reserve(mReads.size() + 1);

  // Add Reference k+1 sequence windows first
  max_num_kmers += (mRegion->Length() - mCurrK + 1);
  labels.emplace_back(Label::REFERENCE);
  kplus_ones.emplace_back(SlidingView(mRegion->SeqView(), mCurrK + 1));

  // Add sample read k+1 sequence windows
  for (const auto& read : mReads) {
    nsample_bases += read.Length();
    max_num_kmers += (read.Length() - mCurrK + 1);
    labels.emplace_back(read.TagKind());
    kplus_ones.emplace_back(SlidingView(read.SeqView(), mCurrK + 1));
  }

  mRefNodeIds.clear();
  mAverageCov = static_cast<f64>(nsample_bases) / static_cast<f64>(mRegion->Length());
  const absl::FixedArray<SeqNodes> added_nodes = AddToGraph(kplus_ones, labels, max_num_kmers);
  const SeqNodes& ref_nodes = added_nodes[0];
  mRefNodeIds.reserve(ref_nodes.size());
  std::transform(ref_nodes.cbegin(), ref_nodes.cend(), std::back_inserter(mRefNodeIds),
                 [](const Node* rnode) -> NodeID { return rnode->Identifier(); });

  mate_mers.clear();

  static constexpr u8 MIN_KMER_BASE_QUALITY = 20;
  static const auto is_low_qual_base = [](const u8 base_qual) -> bool { return base_qual < MIN_KMER_BASE_QUALITY; };

  for (usize rd_idx = 0; rd_idx < mReads.size(); ++rd_idx) {
    const auto read_label = mReads[rd_idx].SrcLabel();
    const auto mm_label = fmt::format("{}{}", mReads[rd_idx].QnameView(), read_label.GetData());
    const auto quals_view = SlidingView(mReads[rd_idx].QualView(), mCurrK);

    for (usize kmer_idx = 0; kmer_idx < added_nodes[rd_idx + 1].size(); ++kmer_idx) {
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (std::ranges::any_of(quals_view[kmer_idx], is_low_qual_base)) continue;

      auto* node = added_nodes[rd_idx + 1][kmer_idx];
      auto mm_pair = std::make_pair(mm_label, node->Identifier());
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (mate_mers.contains(mm_pair)) return;

      node->IncrementReadSupport(read_label);
      mate_mers.emplace(std::move(mm_pair));
    }
  }
}

auto Graph::AddToGraph(SeqKplusOnes kplus_ones, SeqLabels labels, const usize max_kmers) -> GraphNodes {
  mNodes.clear();
  mNodes.reserve(max_kmers);
  absl::FixedArray<SeqNodes> results(kplus_ones.size());
  static constexpr auto dflt_order = Kmer::Ordering::DEFAULT;

  for (usize seq_idx = 0; seq_idx < kplus_ones.size(); ++seq_idx) {
    const auto& curr_label = labels[seq_idx];
    results[seq_idx].reserve(kplus_ones[seq_idx].size() + 1);

    for (usize mer_idx = 0; mer_idx < kplus_ones[seq_idx].size(); ++mer_idx) {
      const auto seq1 = absl::ClippedSubstr(kplus_ones[seq_idx][mer_idx], 0, mCurrK);
      const auto seq2 = absl::ClippedSubstr(kplus_ones[seq_idx][mer_idx], 1, mCurrK);

      const auto left_id = CanonicalKmerHash(seq1);
      const auto right_id = CanonicalKmerHash(seq2);

      mNodes.try_emplace(left_id, std::make_unique<Node>(seq1, curr_label));
      mNodes.try_emplace(right_id, std::make_unique<Node>(seq2, curr_label));

      auto& first = mNodes.at(left_id);
      auto& second = mNodes.at(right_id);

      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (mer_idx == 0) results[seq_idx].emplace_back(first.get());

      const auto fwd_edge = MakeFwdEdgeKind({first->SignFor(dflt_order), second->SignFor(dflt_order)});
      first->EmplaceEdge(NodeIDPair{left_id, right_id}, fwd_edge);
      second->EmplaceEdge(NodeIDPair{right_id, left_id}, RevEdgeKind(fwd_edge));

      results[seq_idx].emplace_back(second.get());
    }
  }

  return results;
}

auto Graph::CanonicalKmerHash(std::string_view seq) -> u64 {
  auto rc_seq = RevComp(seq);
  return seq < rc_seq ? HashStr64(seq) : HashStr64(rc_seq);
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

void Graph::WriteDot(State state, usize comp_id) {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mParams.mOutGraphsDir.empty()) return;

  using namespace std::string_view_literals;
  const auto win_id = fmt::format("{}_{}_{}", mRegion->ChromName(), mRegion->StartPos1(), mRegion->EndPos1());
  const auto fname = fmt::format("dbg__{}__{}__k{}__comp{}.dot", win_id, ToString(state), mCurrK, comp_id);

  const auto out_path = mParams.mOutGraphsDir / "dbg_graph" / fname;
  std::filesystem::create_directories(mParams.mOutGraphsDir / "dbg_graph");
  return SerializeToDot(mNodes, out_path, comp_id, {mSourceAndSinkIds.cbegin(), mSourceAndSinkIds.cend()});
}

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
    const auto fill_color = is_background_node                     ? "darkgray"sv
                            : nodes_highlight.contains(item.first) ? "orchid"sv
                            : item.second->IsShared()              ? "steelblue"sv
                            : item.second->IsTumorOnly()           ? "indianred"sv
                            : item.second->IsNormalOnly()          ? "mediumseagreen"sv
                                                                   : "lightblue"sv;

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
