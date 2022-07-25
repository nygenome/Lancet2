#include "lancet2/node.h"

#include <algorithm>
#include <array>
#include <limits>

#include "lancet2/assert_macro.h"
#include "lancet2/utils.h"

namespace lancet2 {
Node::Node(const Kmer& k) : mer(k), nodeID(k.GetHash()), quals(k.GetLength()), labels(k.GetLength()) {}

auto Node::CanMerge(const Node& buddy, BuddyPosition merge_dir, usize k) const -> bool {
  if (IsMockNode() || buddy.IsMockNode()) return false;
  const auto reverseBuddy = buddy.GetOrientation() != GetOrientation();
  return mer.CanMergeKmers(buddy.mer, merge_dir, reverseBuddy, k);
}

void Node::MergeBuddy(const Node& buddy, BuddyPosition dir, usize k) {
  // Everything except edges are merged from buddy into the node
  const auto reverseBuddy = buddy.GetOrientation() != GetOrientation();
  Reserve(mer.GetLength() + buddy.GetLength() - k + 1);

  mer.MergeBuddy(buddy.mer, dir, reverseBuddy, k);
  quals.MergeBuddy(buddy.quals, dir, reverseBuddy, k);
  covs.MergeBuddy(buddy.covs, mer.GetSize(), buddy.GetLength(), k);
  labels.MergeBuddy(buddy.labels, dir, reverseBuddy, k);

  if (!bxData.IsEmpty() || !buddy.bxData.IsEmpty()) bxData.Merge(buddy.bxData);
}

void Node::EmplaceEdge(NodeIdentifier dest_id, EdgeKind k) {
  if (dest_id == MOCK_SOURCE_ID || dest_id == MOCK_SINK_ID) numMockEdges++;
  if (dest_id == nodeID) numSelfEdges++;

  const auto result = edgeSet.emplace(dest_id, k);
  if (result.second) {
    orderedEdges.emplace_back(dest_id, k);
    std::sort(orderedEdges.begin(), orderedEdges.end());
  }
}

void Node::EraseEdge(NodeIdentifier dest_id, EdgeKind k) {
  auto itr = edgeSet.find(Edge(dest_id, k));
  if (itr != edgeSet.end()) {
    edgeSet.erase(itr);
    orderedEdges.clear();
    orderedEdges.insert(orderedEdges.end(), edgeSet.cbegin(), edgeSet.cend());
    std::sort(orderedEdges.begin(), orderedEdges.end());
  }
}

static constexpr auto ALL_EDGE_KINDS = std::array<EdgeKind, 4>{EdgeKind::FF, EdgeKind::FR, EdgeKind::RF, EdgeKind::RR};

void Node::EraseEdge(NodeIdentifier dest_id) {
  std::for_each(ALL_EDGE_KINDS.cbegin(), ALL_EDGE_KINDS.cend(), [&](const EdgeKind& ek) { EraseEdge(dest_id, ek); });
}

void Node::ClearEdges() {
  orderedEdges.clear();
  edgeSet.clear();
}

auto Node::NumEdges(Strand direction) const -> usize {
  return std::count_if(edgeSet.cbegin(), edgeSet.cend(), [&direction](const Edge& e) {
    // faux nodes are not supported by real data, they exist only for path travesal, skip them in counts
    return e.GetSrcDir() == direction && e.GetDstID() != MOCK_SOURCE_ID && e.GetDstID() != MOCK_SINK_ID;
  });
}

auto Node::NumEdges() const -> usize { return edgeSet.size() - numMockEdges; }

void Node::UpdateQual(std::string_view sv) { quals.Push(sv); }
void Node::UpdateLabel(KmerLabel label) { labels.Push(label); }

void Node::UpdateCovInfo(const ReadInfo& ri, bool is_tenx_mode) {
  if (is_tenx_mode) return covs.Update(BXCount(ri.label, ri.strand), ri.label, ri.strand);
  return covs.Update(ri.label, ri.strand);
}
void Node::IncrementCov(SampleLabel label, Strand s) { covs.Update(label, s); }

auto Node::FillColor() const -> std::string { return IsSource() ? "cyan3" : IsSink() ? "yellow2" : labels.FillColor(); }

auto Node::LabelRatio(KmerLabel label) const -> double { return labels.LabelRatio(label); }
auto Node::HasLabel(KmerLabel label) const -> bool { return labels.HasLabel(label); }
auto Node::IsLabelOnly(KmerLabel label) const -> bool { return labels.IsLabelOnly(label); }

auto Node::TotalSampleCount() const -> u16 {
  return SampleCount(SampleLabel::TUMOR) + SampleCount(SampleLabel::NORMAL);
}
auto Node::SampleCount(SampleLabel label) const -> u16 { return covs.TotalCov(label); }
auto Node::SampleCount(SampleLabel label, Strand s) const -> u16 { return covs.StrandCov(label, s); }
auto Node::BXCount(SampleLabel label, Strand s) const -> u16 { return bxData.BXCount(label, s); }

auto Node::LowQualPositions(u32 min_bq) const -> std::vector<bool> {
  return quals.LowQualPositions(static_cast<double>(min_bq));
}

auto Node::FindMergeableNeighbours() const -> std::vector<NodeNeighbour> {
  if (numSelfEdges != 0 || orderedEdges.size() != 2) return {};

  std::vector<NodeNeighbour> results;
  std::for_each(orderedEdges.cbegin(), orderedEdges.cend(), [&results](const Edge& e) {
    if (e.GetDstID() == MOCK_SOURCE_ID || e.GetDstID() == MOCK_SINK_ID) return;
    results.emplace_back(e);
  });

  return results;
}
}  // namespace lancet2
