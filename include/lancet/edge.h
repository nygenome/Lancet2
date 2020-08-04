#pragma once

#include <cstddef>
#include <cstdint>

#include "lancet/core_enums.h"

namespace lancet {
using EdgeIdentifier = std::size_t;

class Edge {
 public:
  Edge(std::size_t destination_id, EdgeKind kind) : destId(destination_id), edgeKind(kind) {}
  Edge() = default;

  [[nodiscard]] auto ID() const -> std::uint64_t;
  [[nodiscard]] auto Kind() const -> EdgeKind { return edgeKind; }
  [[nodiscard]] auto DestinationID() const -> EdgeIdentifier { return destId; }

  auto operator==(const Edge& other) const -> bool { return destId == other.destId && edgeKind == other.edgeKind; }
  auto operator!=(const Edge& other) const -> bool { return destId != other.destId || edgeKind != other.edgeKind; }

  auto operator<(const Edge& other) const -> bool {
    return static_cast<std::uint8_t>(edgeKind) < static_cast<std::uint8_t>(other.edgeKind) || destId < other.destId;
  }

  [[nodiscard]] auto SrcDirection() const noexcept -> Strand {
    return (edgeKind == EdgeKind::FF || edgeKind == EdgeKind::FR) ? Strand::FWD : Strand::REV;
  }

  [[nodiscard]] auto DstDirection() const noexcept -> Strand {
    return (edgeKind == EdgeKind::FF || edgeKind == EdgeKind::RF) ? Strand::FWD : Strand::REV;
  }

  template <typename H>
  friend auto AbslHashValue(H h, const Edge& e) -> H {
    return H::combine(std::move(h), e.destId, static_cast<std::uint8_t>(e.edgeKind));
  }

 private:
  EdgeIdentifier destId = 0;
  EdgeKind edgeKind = EdgeKind::FF;
};
}  // namespace lancet
