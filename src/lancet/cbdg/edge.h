#ifndef SRC_LANCET_CBDG_EDGE_H_
#define SRC_LANCET_CBDG_EDGE_H_

#include <array>
#include <utility>

#include "lancet/base/types.h"
#include "lancet/cbdg/kmer.h"

namespace lancet::cbdg {

class Edge {
 public:
  Edge() = default;
  explicit Edge(std::array<u64, 2> src_dst_ids, EdgeKind kind)
      : mSrcId(src_dst_ids[0]), mDstId(src_dst_ids[1]), mEdgeKind(kind) {}

  [[nodiscard]] auto SrcId() const noexcept -> u64 { return mSrcId; }
  [[nodiscard]] auto DstId() const noexcept -> u64 { return mDstId; }
  [[nodiscard]] auto Kind() const noexcept -> EdgeKind { return mEdgeKind; }

  [[nodiscard]] auto SrcSign() const noexcept -> Kmer::Sign { return SplitIntoSignPair(mEdgeKind)[0]; }
  [[nodiscard]] auto DstSign() const noexcept -> Kmer::Sign { return SplitIntoSignPair(mEdgeKind)[1]; }

  [[nodiscard]] auto IsSelfLoop() const noexcept -> bool { return mSrcId == mDstId; }
  [[nodiscard]] auto IsSelfMirror() const noexcept -> bool {
    // Any 2 nodes connected to each other have 2 unique edges connecting them. src --> dst and dst --> src
    // The only exception to this rule are self-loops with +- or -+ edge, where there is only one edge.
    return IsSelfLoop() && (mEdgeKind == EdgeKind::PLUS_MINUS || mEdgeKind == EdgeKind::MINUS_PLUS);
  }

  [[nodiscard]] auto MirrorEdge() const noexcept -> Edge { return Edge({mDstId, mSrcId}, RevEdgeKind(mEdgeKind)); }

  template <typename HashState>
  friend auto AbslHashValue(HashState hash_state, const Edge& edge) -> HashState {
    return HashState::combine(std::move(hash_state), edge.mSrcId, edge.mDstId, static_cast<u64>(edge.mEdgeKind));
  }

  auto operator==(const Edge& rhs) const -> bool {
    return mSrcId == rhs.mSrcId && mDstId == rhs.mDstId && mEdgeKind == rhs.mEdgeKind;
  }
  auto operator!=(const Edge& rhs) const -> bool { return !(rhs == *this); }

 private:
  u64 mSrcId = 0;
  u64 mDstId = 0;
  EdgeKind mEdgeKind = EdgeKind::PLUS_PLUS;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_EDGE_H_
