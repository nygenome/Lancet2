#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "lancet/barcode_set.h"
#include "lancet/core_enums.h"
#include "lancet/edge.h"
#include "lancet/kmer.h"
#include "lancet/node_cov.h"
#include "lancet/node_hp.h"
#include "lancet/node_label.h"
#include "lancet/node_neighbour.h"
#include "lancet/node_qual.h"
#include "lancet/read_info.h"

namespace lancet {
using NodeIdentifier = std::size_t;

class Node {
 public:
  std::size_t ComponentID = 0;  // NOLINT

  explicit Node(const Kmer& k);
  explicit Node(NodeIdentifier nid) : nodeID(nid), quals(0), covs(0), labels(0) {}
  Node() = delete;

  [[nodiscard]] auto CanMerge(const Node& buddy, BuddyPosition merge_dir, std::size_t k) const -> bool;

  // Everything except edges are merged from `buddy` into the node
  // NOTE: nodeID is not updated for merged nodes. Original node ID is maintained so the hash table is still valid
  void MergeBuddy(const Node& buddy, BuddyPosition dir, std::size_t k);

  [[nodiscard]] auto FwdSeq() const noexcept -> std::string { return mer.FwdSeq(); }
  [[nodiscard]] auto SeqView() const noexcept -> std::string_view { return mer.SeqView(); }
  [[nodiscard]] auto Orientation() const noexcept -> Strand { return mer.Orientation(); }
  [[nodiscard]] auto Length() const noexcept -> std::size_t { return mer.Length(); }

  [[nodiscard]] auto IsEmpty() const noexcept -> bool {
    return mer.IsEmpty() && quals.IsEmpty() && covs.IsEmpty() && labels.IsEmpty();
  }

  [[nodiscard]] auto ID() const -> NodeIdentifier { return nodeID; }
  [[nodiscard]] auto IsSource() const -> bool { return nodeID == MOCK_SOURCE_ID; }
  [[nodiscard]] auto IsSink() const -> bool { return nodeID == MOCK_SINK_ID; }
  [[nodiscard]] auto IsMockNode() const -> bool { return IsSource() || IsSink(); }

  friend auto operator==(const Node& lhs, const Node& rhs) -> bool { return lhs.mer.ID() == rhs.mer.ID(); }
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

  [[nodiscard]] auto NumEdges(Strand direction) const -> std::size_t;
  [[nodiscard]] auto NumEdges() const -> std::size_t;

  void UpdateQual(std::string_view sv);
  void UpdateLabel(KmerLabel label);

  void UpdateHPInfo(const ReadInfo& ri, std::uint32_t min_base_qual);
  void UpdateCovInfo(const ReadInfo& ri, std::uint32_t min_base_qual, bool is_tenx_mode);
  void IncrementCov(SampleLabel label, Strand s, std::size_t base_position);

  [[nodiscard]] auto FillColor() const -> std::string;

  [[nodiscard]] auto LabelRatio(KmerLabel label) const -> double;
  [[nodiscard]] auto HasLabel(KmerLabel label) const -> bool;
  [[nodiscard]] auto IsLabelOnly(KmerLabel label) const -> bool;

  [[nodiscard]] auto TotalSampleCount() const -> std::uint16_t;
  [[nodiscard]] auto SampleCount(SampleLabel label) const -> std::uint16_t;
  [[nodiscard]] auto SampleCount(SampleLabel label, Strand s) const -> std::uint16_t;
  [[nodiscard]] auto BXCount(SampleLabel label, Strand s) const -> std::uint16_t;

  [[nodiscard]] auto CovAt(SampleLabel label, std::size_t pos) -> BaseCov& { return covs.At(label, pos); }
  [[nodiscard]] auto CovAt(SampleLabel label, std::size_t pos) const -> const BaseCov& { return covs.At(label, pos); }
  [[nodiscard]] auto MinSampleBaseCov(bool bq_pass = false) const -> std::uint16_t;

  [[nodiscard]] auto LowQualPositions(std::uint32_t min_bq) const -> std::vector<bool>;

  [[nodiscard]] auto CovData() const noexcept -> NodeCov { return covs; }
  [[nodiscard]] auto HasBXData() const noexcept -> bool { return !bxData.IsEmpty(); }
  [[nodiscard]] auto hasHPData() const noexcept -> bool { return !hpData.IsEmpty(); }

  [[nodiscard]] auto HPData() const noexcept -> NodeHP { return hpData.IsEmpty() ? NodeHP(covs) : hpData; }

  [[nodiscard]] auto FindMergeableNeighbours() const -> std::vector<NodeNeighbour>;

  void Reserve(const std::size_t count) {
    quals.Reserve(count);
    covs.Reserve(count);
    labels.Reserve(count);
    hpData.Reserve(count);
  }

 private:
  Kmer mer;
  NodeIdentifier nodeID;

  std::size_t numMockEdges = 0;
  std::size_t numSelfEdges = 0;
  EdgeContainer orderedEdges;
  absl::flat_hash_set<Edge> edgeSet;

  NodeQual quals;
  NodeCov covs;
  NodeLabel labels;

  BarcodeSet bxData;
  NodeHP hpData;
};  // namespace lancet
}  // namespace lancet
