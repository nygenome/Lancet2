#pragma once

#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "lancet2/barcode_set.h"
#include "lancet2/core_enums.h"
#include "lancet2/edge.h"
#include "lancet2/kmer.h"
#include "lancet2/node_cov.h"
#include "lancet2/node_hp.h"
#include "lancet2/node_label.h"
#include "lancet2/node_neighbour.h"
#include "lancet2/node_qual.h"
#include "lancet2/read_info.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
using NodeIdentifier = usize;

class Node {
 public:
  usize ComponentID = 0;  // NOLINT

  explicit Node(const Kmer& k);
  explicit Node(NodeIdentifier nid) : nodeID(nid), quals(0), covs(0), labels(0) {}
  Node() = delete;

  [[nodiscard]] auto CanMerge(const Node& buddy, BuddyPosition merge_dir, usize k) const -> bool;

  // Everything except edges are merged from `buddy` into the node
  // NOTE: nodeID is not updated for merged nodes. Original node ID is maintained so the hash table is still valid
  void MergeBuddy(const Node& buddy, BuddyPosition dir, usize k);

  [[nodiscard]] auto GetFwdSeq() const noexcept -> std::string { return mer.GetFwdSeq(); }
  [[nodiscard]] auto GetSeqView() const noexcept -> std::string_view { return mer.GetSeqView(); }
  [[nodiscard]] auto GetOrientation() const noexcept -> Strand { return mer.GetOrientation(); }
  [[nodiscard]] auto GetLength() const noexcept -> usize { return mer.GetLength(); }

  [[nodiscard]] auto IsEmpty() const noexcept -> bool {
    return mer.IsEmpty() && quals.IsEmpty() && covs.IsEmpty() && labels.IsEmpty();
  }

  [[nodiscard]] auto GetID() const -> NodeIdentifier { return nodeID; }
  [[nodiscard]] auto IsSource() const -> bool { return nodeID == MOCK_SOURCE_ID; }
  [[nodiscard]] auto IsSink() const -> bool { return nodeID == MOCK_SINK_ID; }
  [[nodiscard]] auto IsMockNode() const -> bool { return IsSource() || IsSink(); }

  friend auto operator==(const Node& lhs, const Node& rhs) -> bool { return lhs.mer.GetHash() == rhs.mer.GetHash(); }
  friend auto operator!=(const Node& lhs, const Node& rhs) -> bool { return !(lhs == rhs); }

  using EdgeContainer = std::vector<Edge>;
  using EdgeIterator = EdgeContainer ::iterator;
  using ConstEdgeIterator = EdgeContainer::const_iterator;

  [[nodiscard]] auto begin() -> EdgeIterator { return orderedEdges.begin(); }
  [[nodiscard]] auto end() -> EdgeIterator { return orderedEdges.end(); }
  [[nodiscard]] auto begin() const -> ConstEdgeIterator { return orderedEdges.begin(); }
  [[nodiscard]] auto cbegin() const -> ConstEdgeIterator { return orderedEdges.cbegin(); }
  [[nodiscard]] auto end() const -> ConstEdgeIterator { return orderedEdges.end(); }
  [[nodiscard]] auto cend() const -> ConstEdgeIterator { return orderedEdges.cend(); }

  void EmplaceEdge(NodeIdentifier dest_id, EdgeKind k);
  void EraseEdge(NodeIdentifier dest_id, EdgeKind k);
  void EraseEdge(NodeIdentifier dest_id);

  void ClearEdges();
  [[nodiscard]] auto HasSelfLoop() const -> bool;
  [[nodiscard]] auto HasConnection(NodeIdentifier dest) const -> bool;

  [[nodiscard]] auto NumEdges(Strand direction) const -> usize;
  [[nodiscard]] auto NumEdges() const -> usize;

  void UpdateQual(std::string_view sv);
  void UpdateLabel(KmerLabel label);

  void UpdateHPInfo(const ReadInfo& ri, u32 min_base_qual);
  void UpdateCovInfo(const ReadInfo& ri, u32 min_base_qual, bool is_tenx_mode);
  void IncrementCov(SampleLabel label, Strand s, usize base_position);

  [[nodiscard]] auto FillColor() const -> std::string;

  [[nodiscard]] auto TumorOnlyLabelRatio(KmerLabel label) const -> double;
  [[nodiscard]] auto HasLabel(KmerLabel label) const -> bool;
  [[nodiscard]] auto IsLabelOnly(KmerLabel label) const -> bool;

  [[nodiscard]] auto TotalSampleCount() const -> u16;
  [[nodiscard]] auto SampleCount(SampleLabel label) const -> u16;
  [[nodiscard]] auto SampleCount(SampleLabel label, Strand s) const -> u16;
  [[nodiscard]] auto BXCount(SampleLabel label, Strand s) const -> u16;

  [[nodiscard]] auto CovAt(SampleLabel label, usize pos) -> BaseCov& { return covs.At(label, pos); }
  [[nodiscard]] auto CovAt(SampleLabel label, usize pos) const -> const BaseCov& { return covs.At(label, pos); }
  [[nodiscard]] auto MinSampleBaseCov(bool bq_pass = false) const -> u16;

  [[nodiscard]] auto LowQualPositions(u32 min_bq) const -> std::vector<bool>;

  [[nodiscard]] auto CovData() const noexcept -> NodeCov { return covs; }
  [[nodiscard]] auto HasBXData() const noexcept -> bool { return !bxData.IsEmpty(); }
  [[nodiscard]] auto hasHPData() const noexcept -> bool { return !hpData.IsEmpty(); }

  [[nodiscard]] auto HPData() const noexcept -> NodeHP { return hpData.IsEmpty() ? NodeHP(covs) : hpData; }

  [[nodiscard]] auto FindMergeableNeighbours() const -> std::vector<NodeNeighbour>;

  void Reserve(const usize count) {
    quals.Reserve(count);
    covs.Reserve(count);
    labels.Reserve(count);
    hpData.Reserve(count);
  }

 private:
  Kmer mer;
  NodeIdentifier nodeID;

  usize numMockEdges = 0;
  usize numSelfEdges = 0;
  EdgeContainer orderedEdges;
  absl::flat_hash_set<Edge> edgeSet;

  NodeQual quals;
  NodeCov covs;
  NodeLabel labels;

  BarcodeSet bxData;
  NodeHP hpData;
};  // namespace lancet2
}  // namespace lancet2
