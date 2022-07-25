#pragma once

#include <string>
#include <string_view>

#include "absl/container/fixed_array.h"
#include "absl/types/span.h"
#include "lancet2/core_enums.h"
#include "lancet2/node.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
struct EdgeNodeIds {
  NodeIdentifier srcId = 0;
  NodeIdentifier dstId = 0;
};

using PathNodeIds = absl::FixedArray<EdgeNodeIds>;

class Path {
 public:
  Path(absl::FixedArray<const Node*> path_nodes, absl::FixedArray<const Edge*> path_edges, std::string path_seq);
  Path() = delete;

  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return pathSeq.empty(); }
  [[nodiscard]] auto GetLength() const noexcept -> usize { return pathSeq.length(); }
  [[nodiscard]] auto GetNumNodes() const noexcept -> usize { return nodesList.size(); }
  [[nodiscard]] auto GetSeqView() const noexcept -> std::string_view { return pathSeq; }

  [[nodiscard]] auto FindSpanningNode(usize path_pos, usize curr_k) const -> const Node*;

  auto operator==(const Path& other) const -> bool { return pathSeq == other.pathSeq; }
  auto operator!=(const Path& other) const -> bool { return !(*this == other); }

  [[nodiscard]] auto TouchedEdges() const -> absl::Span<const Edge* const> { return absl::MakeConstSpan(edgesList); }
  [[nodiscard]] auto TouchedNodes() const -> absl::Span<const Node* const> { return absl::MakeConstSpan(nodesList); }

  [[nodiscard]] auto TouchedEdgeIDs() const -> PathNodeIds;

 private:
  absl::FixedArray<const Node*> nodesList;
  absl::FixedArray<const Edge*> edgesList;
  std::string pathSeq;
};
}  // namespace lancet2
