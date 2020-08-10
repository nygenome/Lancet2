#pragma once

#include <cstddef>
#include <string>
#include <string_view>

#include "absl/container/fixed_array.h"
#include "absl/types/span.h"
#include "lancet/base_hpcov.h"
#include "lancet/core_enums.h"
#include "lancet/node.h"

namespace lancet {
struct EdgeNodeIds {
  NodeIdentifier srcId = 0;
  NodeIdentifier dstId = 0;
};

using PathNodeIds = absl::FixedArray<EdgeNodeIds>;

class Path {
 public:
  Path(absl::FixedArray<const Node*> path_nodes, absl::FixedArray<const Edge*> path_edges, std::string path_seq,
       NodeCov path_covs, NodeHP path_hps);
  Path() = delete;

  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return pathSeq.empty(); }
  [[nodiscard]] auto Length() const noexcept -> std::size_t { return pathSeq.length(); }
  [[nodiscard]] auto NumNodes() const noexcept -> std::size_t { return nodesList.size(); }
  [[nodiscard]] auto SeqView() const noexcept -> std::string_view { return pathSeq; }

  [[nodiscard]] auto FindSpanningNode(std::size_t path_pos, std::size_t curr_k) const -> const Node*;

  [[nodiscard]] auto HpCovAt(SampleLabel label, std::size_t pos) const -> BaseHpCov {
    if (pathHPs.IsEmpty()) return BaseHpCov{pathCovs.At(label, pos)};
    return BaseHpCov{pathCovs.At(label, pos), pathHPs.At(label, pos)};
  }

  auto operator==(const Path& other) const -> bool { return pathSeq == other.pathSeq; }
  auto operator!=(const Path& other) const -> bool { return !(*this == other); }

  [[nodiscard]] auto TouchedEdges() const -> absl::Span<const Edge* const> { return absl::MakeConstSpan(edgesList); }
  [[nodiscard]] auto TouchedNodes() const -> absl::Span<const Node* const> { return absl::MakeConstSpan(nodesList); }

  [[nodiscard]] auto TouchedEdgeIDs() const -> PathNodeIds;

 private:
  absl::FixedArray<const Node*> nodesList;
  absl::FixedArray<const Edge*> edgesList;
  std::string pathSeq;
  NodeCov pathCovs;
  NodeHP pathHPs;
};
}  // namespace lancet
