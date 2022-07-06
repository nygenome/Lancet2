#pragma once

#include <cstddef>

#include "lancet2/edge.h"

namespace lancet2 {
using NodeIdentifier = std::size_t;

struct NodeNeighbour {
  explicit NodeNeighbour(const Edge& e) : buddyId(e.DestinationID()), edgeKind(e.Kind()) {}
  explicit NodeNeighbour(NodeIdentifier nid, EdgeKind ek) : buddyId(nid), edgeKind(ek) {}

  NodeIdentifier buddyId = 0;        // NOLINT
  EdgeKind edgeKind = EdgeKind::FF;  // NOLINT

  auto operator==(const NodeNeighbour& other) const -> bool {
    return buddyId == other.buddyId && edgeKind == other.edgeKind;
  }
  auto operator!=(const NodeNeighbour& other) const -> bool { return !(*this == other); }

  auto operator<(const NodeNeighbour& other) const -> bool {
    return static_cast<std::uint8_t>(edgeKind) < static_cast<std::uint8_t>(other.edgeKind) || buddyId < other.buddyId;
  }

  template <typename H>
  friend auto AbslHashValue(H h, const NodeNeighbour& n) -> H {
    return H::combine(std::move(h), n.buddyId, static_cast<std::uint8_t>(n.edgeKind));
  }
};
}  // namespace lancet2
