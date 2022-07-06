#pragma once

#include "lancet2/core_enums.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
using EdgeIdentifier = usize;

class Edge {
 public:
  Edge(usize destination_id, EdgeKind kind) : destId(destination_id), edgeKind(kind) {}
  Edge() = default;

  [[nodiscard]] auto GetID() const -> u64;
  [[nodiscard]] auto GetEdgeKind() const -> EdgeKind { return edgeKind; }
  [[nodiscard]] auto GetDstID() const -> EdgeIdentifier { return destId; }

  auto operator==(const Edge& other) const -> bool { return destId == other.destId && edgeKind == other.edgeKind; }
  auto operator!=(const Edge& other) const -> bool { return destId != other.destId || edgeKind != other.edgeKind; }

  auto operator<(const Edge& other) const -> bool {
    return static_cast<u8>(edgeKind) < static_cast<u8>(other.edgeKind) || destId < other.destId;
  }

  [[nodiscard]] auto GetSrcDir() const noexcept -> Strand {
    return (edgeKind == EdgeKind::FF || edgeKind == EdgeKind::FR) ? Strand::FWD : Strand::REV;
  }

  [[nodiscard]] auto GetDstDir() const noexcept -> Strand {
    return (edgeKind == EdgeKind::FF || edgeKind == EdgeKind::RF) ? Strand::FWD : Strand::REV;
  }

  template <typename H>
  friend auto AbslHashValue(H h, const Edge& e) -> H {
    return H::combine(std::move(h), e.destId, static_cast<u8>(e.edgeKind));
  }

 private:
  EdgeIdentifier destId = 0;
  EdgeKind edgeKind = EdgeKind::FF;
};
}  // namespace lancet2
