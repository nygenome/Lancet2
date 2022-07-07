#pragma once

#include "lancet2/edge.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
using NodeIdentifier = usize;

struct NodeNeighbour {
  explicit NodeNeighbour(const Edge& e) : buddyId(e.GetDstID()), edgeKind(e.GetEdgeKind()) {}
  explicit NodeNeighbour(NodeIdentifier nid, EdgeKind ek) : buddyId(nid), edgeKind(ek) {}

  NodeIdentifier buddyId = 0;        // NOLINT
  EdgeKind edgeKind = EdgeKind::FF;  // NOLINT

  auto operator==(const NodeNeighbour& other) const -> bool {
    return buddyId == other.buddyId && edgeKind == other.edgeKind;
  }
  auto operator!=(const NodeNeighbour& other) const -> bool { return !(*this == other); }

  auto operator<(const NodeNeighbour& other) const -> bool {
    return static_cast<u8>(edgeKind) < static_cast<u8>(other.edgeKind) && buddyId < other.buddyId;
  }

  template <typename H>
  friend auto AbslHashValue(H h, const NodeNeighbour& n) -> H {
    return H::combine(std::move(h), n.buddyId, static_cast<u8>(n.edgeKind));
  }
};
}  // namespace lancet2
