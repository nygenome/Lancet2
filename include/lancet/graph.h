#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <ostream>
#include <string_view>
#include <vector>

#include "absl/container/btree_set.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "lancet/base_hpcov.h"
#include "lancet/cli_params.h"
#include "lancet/core_enums.h"
#include "lancet/node.h"
#include "lancet/path.h"
#include "lancet/ref_window.h"
#include "lancet/transcript.h"
#include "lancet/variant_store.h"

namespace lancet {
class Graph {
 public:
  using NodePtr = std::unique_ptr<Node>;
  using NodeContainer = absl::flat_hash_map<NodeIdentifier, NodePtr>;
  using NodeIterator = NodeContainer::iterator;
  using ConstNodeIterator = NodeContainer::const_iterator;

  Graph(std::shared_ptr<const RefWindow> w, NodeContainer&& data, double avg_cov, std::size_t k,
        std::shared_ptr<const CliParams> p);
  Graph() = delete;

  // 0 = NORMAL, 1 = TUMOR
  using RefInfos = std::array<absl::Span<const BaseHpCov>, 2>;
  void ProcessGraph(RefInfos&& ref_infos, std::vector<Variant>* results);

  [[nodiscard]] auto ShouldIncrementK() const noexcept -> bool { return shouldIncrementK; }

  struct ComponentInfo {
    std::size_t ID = 0;
    std::size_t numNodes = 0;
  };

  [[nodiscard]] auto MarkConnectedComponents() -> std::vector<ComponentInfo>;

  struct SrcSnkResult {
    bool foundSrcAndSnk = false;
    std::size_t startOffset = 0;
    std::size_t endOffset = 0;
  };

  auto MarkSourceSink(std::size_t comp_id) -> SrcSnkResult;
  static auto RefAnchorLen(const SrcSnkResult& r) noexcept -> std::size_t { return r.endOffset - r.startOffset; }

  auto RemoveLowCovNodes(std::size_t comp_id) -> bool;
  auto CompressGraph(std::size_t comp_id) -> bool;
  auto RemoveTips(std::size_t comp_id) -> bool;
  auto RemoveShortLinks(std::size_t comp_id) -> bool;

  [[nodiscard]] auto HasCycle() const -> bool;

  void ProcessPath(const Path& path, const RefInfos& ref_infos, const SrcSnkResult& einfo,
                   std::vector<Variant>* results) const;

  void WriteDot(std::size_t comp_id, const std::string& suffix) const;
  void WriteDot(std::size_t comp_id, absl::Span<const PathNodeIds> flow_paths) const;

  void EraseNode(NodeIterator itr);
  void EraseNode(NodeIdentifier node_id);

  template <class NodeIDIterator>
  void RemoveNodes(NodeIDIterator first, NodeIDIterator last) {
    for (; first != last; ++first) EraseNode(*first);
  }

  void Clear() { return nodesMap.clear(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return nodesMap.empty(); }
  [[nodiscard]] auto Size() const noexcept -> std::size_t { return nodesMap.size(); }

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
  std::size_t kmerSize = 0;
  std::shared_ptr<const CliParams> params = nullptr;
  bool shouldIncrementK = false;
  NodeContainer nodesMap;

  struct RefEndResult {
    NodeIdentifier nodeId = 0;
    std::size_t refMerIdx = 0;
    bool foundEnd = false;
  };

  [[nodiscard]] auto FindRefEnd(GraphEnd k, std::size_t comp_id, absl::Span<const NodeIdentifier> ref_mer_hashes) const
      -> RefEndResult;

  [[nodiscard]] auto FindCompressibleNeighbours(NodeIdentifier src_id) const -> absl::btree_set<NodeNeighbour>;

  void CompressNode(NodeIdentifier src_id, const absl::btree_set<NodeNeighbour>& buddies,
                    absl::flat_hash_set<NodeIdentifier>* compressed) const;

  auto HasCycle(NodeIdentifier node_id, Strand direction, absl::flat_hash_set<NodeIdentifier>* touched) const -> bool;

  static auto ClampToSourceSink(const RefInfos& refs, const SrcSnkResult& ends) -> RefInfos;
  static void ResetSourceSink(const NodeContainer& nc, std::size_t current_component);
  static void DisconnectEdges(NodeIterator itr, const NodeContainer& nc, Strand direction);
};
}  // namespace lancet
