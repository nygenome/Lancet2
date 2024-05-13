#ifndef SRC_LANCET_CBDG_NODE_H_
#define SRC_LANCET_CBDG_NODE_H_

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/label.h"

namespace lancet::cbdg {

using NodeID = u64;
using NodeIDPair = std::array<NodeID, 2>;

class Node {
 public:
  using EdgeSet = absl::flat_hash_set<Edge>;

  Node() = default;
  Node(Kmer&& mer, Label label) : mKmer(std::move(mer)), mLabel(label) {}

  void AddLabel(const Label& label);
  void IncrementReadSupport(const Label& label);

  template <class... Args>
  void EmplaceEdge(Args&&... args) {
    mEdges.emplace(std::forward<Args>(args)...);
  }

  void EraseEdge(const Edge& edge) { mEdges.erase(edge); }
  void EraseAllEdges() { mEdges.clear(); }

  [[nodiscard]] auto NumOutEdges() const noexcept -> usize { return mEdges.size(); }
  [[nodiscard]] auto SeqLength() const noexcept -> usize { return mKmer.Length(); }

  void SetComponentId(const usize comp_id) { mCompId = comp_id; }
  [[nodiscard]] auto GetComponentId() const noexcept -> usize { return mCompId; }

  [[nodiscard]] auto HasTag(const Label::Tag tag) const noexcept -> bool { return mLabel.HasTag(tag); }

  [[nodiscard]] auto IsShared() const noexcept -> bool {
    return HasTag(Label::NORMAL) && HasTag(Label::TUMOR) && !HasTag(Label::REFERENCE);
  }

  [[nodiscard]] auto IsNormalOnly() const noexcept -> bool {
    return HasTag(Label::NORMAL) && !HasTag(Label::TUMOR) && !HasTag(Label::REFERENCE);
  }

  [[nodiscard]] auto IsTumorOnly() const noexcept -> bool {
    return HasTag(Label::TUMOR) && !HasTag(Label::NORMAL) && !HasTag(Label::REFERENCE);
  }

  [[nodiscard]] auto NormalReadSupport() const noexcept -> u32;
  [[nodiscard]] auto TumorReadSupport() const noexcept -> u32;
  [[nodiscard]] auto TotalReadSupport() const noexcept -> u32;

  [[nodiscard]] auto KmerData() const noexcept -> Kmer { return mKmer; }
  [[nodiscard]] auto Identifier() const noexcept -> NodeID { return mKmer.Identifier(); }
  [[nodiscard]] auto Length() const noexcept -> usize { return mKmer.Length(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return mKmer.IsEmpty(); }

  [[nodiscard]] auto SignFor(const Kmer::Ordering ord) const noexcept -> Kmer::Sign { return mKmer.SignFor(ord); }
  [[nodiscard]] auto SequenceFor(const Kmer::Ordering ord) const -> std::string { return mKmer.SequenceFor(ord); }

  void Merge(const Node& other, EdgeKind conn_kind, usize currk);

  [[nodiscard]] auto HasSelfLoop() const -> bool;
  [[nodiscard]] auto FindEdgesInDirection(Kmer::Ordering ord) const -> std::vector<Edge>;

  using EdgeIterator = EdgeSet::iterator;
  using EdgeConstIterator = EdgeSet::const_iterator;

  [[nodiscard]] auto begin() -> EdgeIterator { return mEdges.begin(); }
  [[nodiscard]] auto end() -> EdgeIterator { return mEdges.end(); }

  [[nodiscard]] auto begin() const -> EdgeConstIterator { return mEdges.begin(); }
  [[nodiscard]] auto end() const -> EdgeConstIterator { return mEdges.end(); }
  [[nodiscard]] auto cbegin() const -> EdgeConstIterator { return mEdges.cbegin(); }
  [[nodiscard]] auto cend() const -> EdgeConstIterator { return mEdges.cend(); }

 private:
  static constexpr usize NORMAL_COUNT_INDEX = 0;
  static constexpr usize TUMOR_COUNT_INDEX = 1;
  using Counts = std::array<u32, 2>;

  Kmer mKmer;
  Label mLabel;
  EdgeSet mEdges;

  usize mCompId = 0;
  Counts mCounts = {0, 0};
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_NODE_H_
