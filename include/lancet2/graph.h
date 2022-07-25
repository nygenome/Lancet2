#pragma once

#include <memory>
#include <ostream>
#include <string_view>
#include <vector>

#include "absl/container/btree_set.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "lancet2/base_hpcov.h"
#include "lancet2/cli_params.h"
#include "lancet2/core_enums.h"
#include "lancet2/node.h"
#include "lancet2/path.h"
#include "lancet2/read_info.h"
#include "lancet2/ref_window.h"
#include "lancet2/sized_ints.h"
#include "lancet2/transcript.h"
#include "lancet2/variant_store.h"

namespace lancet2 {
class Graph {
 public:
  using NodePtr = std::unique_ptr<Node>;
  using NodeContainer = absl::flat_hash_map<NodeIdentifier, NodePtr>;
  using NodeIterator = NodeContainer::iterator;
  using ConstNodeIterator = NodeContainer::const_iterator;

  Graph(std::shared_ptr<const RefWindow> w, NodeContainer&& data, double avg_cov, usize k,
        std::shared_ptr<const CliParams> p);
  Graph() = delete;

  void ProcessGraph(absl::Span<const ReadInfo> reads, std::vector<Variant>* results);

  [[nodiscard]] auto ShouldIncrementK() const noexcept -> bool { return shouldIncrementK; }

  struct ComponentInfo {
    usize ID = 0;
    usize numNodes = 0;
  };

  [[nodiscard]] auto MarkConnectedComponents() -> std::vector<ComponentInfo>;

  struct SrcSnkResult {
    bool foundSrcAndSnk = false;
    usize startOffset = 0;
    usize endOffset = 0;
  };

  auto MarkSourceSink(usize comp_id) -> SrcSnkResult;
  static auto RefAnchorLen(const SrcSnkResult& r) noexcept -> usize { return r.endOffset - r.startOffset; }

  auto RemoveLowCovNodes(usize comp_id) -> bool;
  auto CompressGraph(usize comp_id) -> bool;
  auto RemoveTips(usize comp_id) -> bool;
  auto RemoveShortLinks(usize comp_id) -> bool;

  [[nodiscard]] auto HasCycle() const -> bool;

  void ProcessPath(const Path& path, absl::Span<const ReadInfo> reads, const SrcSnkResult& einfo,
                   std::vector<Variant>* results) const;

  void WritePathFasta(std::string_view path_seq, usize comp_id, usize path_num) const;

  void WriteDot(usize comp_id, const std::string& suffix) const;
  void WriteDot(usize comp_id, absl::Span<const PathNodeIds> flow_paths) const;

  void EraseNode(NodeIterator itr);
  void EraseNode(NodeIdentifier node_id);

  template <class NodeIDIterator>
  void RemoveNodes(NodeIDIterator first, NodeIDIterator last) {
    for (; first != last; ++first) EraseNode(*first);
  }

  void Clear() { return nodesMap.clear(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return nodesMap.empty(); }
  [[nodiscard]] auto Size() const noexcept -> usize { return nodesMap.size(); }

  [[nodiscard]] auto Contains(NodeIdentifier node_id) const -> bool { return nodesMap.contains(node_id); }

  auto begin() -> NodeIterator { return nodesMap.begin(); }
  auto end() -> NodeIterator { return nodesMap.end(); }

  auto begin() const -> ConstNodeIterator { return nodesMap.begin(); }    // NOLINT
  auto cbegin() const -> ConstNodeIterator { return nodesMap.cbegin(); }  // NOLINT
  auto end() const -> ConstNodeIterator { return nodesMap.end(); }        // NOLINT
  auto cend() const -> ConstNodeIterator { return nodesMap.cend(); }      // NOLINT

  [[nodiscard]] auto find(NodeIdentifier node_id) -> NodeIterator { return nodesMap.find(node_id); }
  [[nodiscard]] auto find(NodeIdentifier node_id) const -> ConstNodeIterator { return nodesMap.find(node_id); }

  class DotSerializer;

 private:
  std::shared_ptr<const RefWindow> window;
  double avgSampleCov = 0.0;
  usize kmerSize = 0;
  std::shared_ptr<const CliParams> params = nullptr;
  bool shouldIncrementK = false;
  NodeContainer nodesMap;

  struct RefEndResult {
    NodeIdentifier nodeId = 0;
    usize refMerIdx = 0;
    bool foundEnd = false;
  };

  [[nodiscard]] auto FindRefEnd(GraphEnd k, usize comp_id, absl::Span<const NodeIdentifier> ref_mer_hashes) const
      -> RefEndResult;

  [[nodiscard]] auto FindCompressibleNeighbours(NodeIdentifier src_id) const -> absl::btree_set<NodeNeighbour>;

  void CompressNode(NodeIdentifier src_id, const absl::btree_set<NodeNeighbour>& buddies,
                    absl::flat_hash_set<NodeIdentifier>* compressed) const;

  auto HasCycle(NodeIdentifier node_id, Strand direction, absl::flat_hash_set<NodeIdentifier>* touched) const -> bool;

  static void DisconnectEdgesTo(NodeIterator itr, const NodeContainer& nc);
};
}  // namespace lancet2
