#ifndef SRC_LANCET_CBDG_NODE_H_
#define SRC_LANCET_CBDG_NODE_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/label.h"

#include "absl/container/inlined_vector.h"

#include <algorithm>
#include <array>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lancet::cbdg {

/// Map Label::Tag to a 0-based index into the 2-element mRoleCounts array.
/// CTRL → 0, CASE → 1. Only two roles exist regardless of sample count.
/// REFERENCE is never used for read support tracking.
constexpr auto RoleIndex(Label::Tag const tag) -> usize {
  return tag == Label::CASE ? 1 : 0;
}

using NodeID = u64;
using NodeIDPair = std::array<NodeID, 2>;

// Node stores a canonical k-mer in the colored bidirected de Bruijn graph.
//
// Two tracking systems coexist:
//   Label (3-bit bitmask)     — which roles contributed reads (CTRL|CASE|REF).
//                                Used for graph pruning, DOT visualization, and
//                                somatic state classification.
//   mCounts (per-sample u32)  — how many reads each individual sample contributed.
//                                Used for coverage thresholds and ML features.
//   mRoleCounts (2-element)   — per-role totals (sum of all CTRL samples, sum of
//                                all CASE samples). Fast-path for pruning.
class Node {
 public:
  using EdgeList = absl::InlinedVector<Edge, 8>;

  Node() = default;
  Node(Kmer&& mer, Label label) : mKmer(std::move(mer)), mLabel(label) {}

  void AddLabel(Label const& label);

  /// Increment the per-sample counter and the per-role aggregate.
  /// Grows the per-sample vector on demand if sample_index exceeds current size.
  void IncrementReadSupport(usize sample_index, Label::Tag tag);

  /// Read support for a specific sample. Returns 0 for untracked indices.
  [[nodiscard]] auto ReadSupportForSample(usize sample_index) const -> u32;

  /// Total read support across all samples assigned to the given role.
  [[nodiscard]] auto ReadSupportForRole(Label::Tag role) const noexcept -> u32;

  template <class... Args>
  void EmplaceEdge(Args&&... args) {
    Edge const new_edge(std::forward<Args>(args)...);
    if (std::ranges::find(mEdges, new_edge) == mEdges.end()) {
      mEdges.emplace_back(new_edge);
    }
  }

  void EraseEdge(Edge const& edge) {
    auto* iter = std::ranges::find(mEdges, edge);
    if (iter != mEdges.end()) {
      mEdges.erase(iter);
    }
  }
  void EraseAllEdges() { mEdges.clear(); }

  [[nodiscard]] auto NumOutEdges() const noexcept -> usize { return mEdges.size(); }
  [[nodiscard]] auto SeqLength() const noexcept -> usize { return mKmer.Length(); }

  void SetComponentId(usize const comp_id) { mCompId = comp_id; }
  [[nodiscard]] auto GetComponentId() const noexcept -> usize { return mCompId; }

  [[nodiscard]] auto HasTag(Label::Tag const tag) const noexcept -> bool {
    return mLabel.HasTag(tag);
  }

  [[nodiscard]] auto IsShared() const noexcept -> bool {
    return HasTag(Label::CTRL) && HasTag(Label::CASE) && !HasTag(Label::REFERENCE);
  }

  [[nodiscard]] auto IsCtrlOnly() const noexcept -> bool {
    return HasTag(Label::CTRL) && !HasTag(Label::CASE) && !HasTag(Label::REFERENCE);
  }

  [[nodiscard]] auto IsCaseOnly() const noexcept -> bool {
    return HasTag(Label::CASE) && !HasTag(Label::CTRL) && !HasTag(Label::REFERENCE);
  }

  [[nodiscard]] auto TotalReadSupport() const noexcept -> u32;

  /// Coverage-bounded confidence score for walk enumeration and SPOA weighting.
  ///
  /// confidence = floor(TotalReadSupport × concordance) + (REFERENCE ? 1 : 0)
  ///   concordance = confirming_samples / num_samples  ∈ (0, 1]
  ///
  /// @param num_samples Authoritative total sample count from GraphParams.
  ///   Cannot be derived from mCounts.size() because IncrementReadSupport
  ///   grows the vector lazily — nodes only touched by low-index samples
  ///   have an undersized mCounts, which would inflate concordance to 1.0.
  ///
  /// Result ∈ [0, TotalReadSupport + 1]:
  ///   - Maximum (TotalReadSupport + 1) when all samples confirm AND k-mer is in the reference.
  ///   - Zero when TotalReadSupport is zero.
  ///   - One for singleton nodes (every sample has ≤1 read).
  ///
  /// The +1 reference bonus is additive, not multiplicative — it gives
  /// reference-confirmed k-mers a tiebreaker without distorting the coverage
  /// scale. This keeps REF and ALT path weights directly comparable.
  [[nodiscard]] auto Confidence(usize num_samples) const noexcept -> u32;

  /// True if every sample with non-zero coverage has exactly 1 read.
  /// N-sample generalization of the 2-sample (ctrl==1 && case==1) check.
  [[nodiscard]] auto IsAllSingletons() const noexcept -> bool;

  [[nodiscard]] auto KmerData() const noexcept -> Kmer const& { return mKmer; }
  [[nodiscard]] auto Identifier() const noexcept -> NodeID { return mKmer.Identifier(); }
  [[nodiscard]] auto Length() const noexcept -> usize { return mKmer.Length(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return mKmer.IsEmpty(); }

  [[nodiscard]] auto SignFor(Kmer::Ordering const ord) const noexcept -> Kmer::Sign {
    return mKmer.SignFor(ord);
  }
  [[nodiscard]] auto SequenceFor(Kmer::Ordering const ord) const -> std::string {
    return mKmer.SequenceFor(ord);
  }

  /// Non-owning view of the canonical k-mer sequence (zero allocation).
  /// Forwards to `Kmer::SeqView()`. Use on hot paths instead of
  /// `SequenceFor(DEFAULT)` which allocates a fresh string per call.
  [[nodiscard]] auto SeqView() const noexcept -> std::string_view { return mKmer.SeqView(); }

  void Merge(Node const& other, EdgeKind conn_kind, usize currk);

  [[nodiscard]] auto HasSelfLoop() const -> bool;
  [[nodiscard]] auto FindEdgesInDirection(Kmer::Ordering ord) const -> std::vector<Edge>;

  using EdgeIterator = EdgeList::iterator;
  using EdgeConstIterator = EdgeList::const_iterator;

  [[nodiscard]] auto begin() -> EdgeIterator { return mEdges.begin(); }
  [[nodiscard]] auto end() -> EdgeIterator { return mEdges.end(); }

  [[nodiscard]] auto begin() const -> EdgeConstIterator { return mEdges.begin(); }
  [[nodiscard]] auto end() const -> EdgeConstIterator { return mEdges.end(); }
  [[nodiscard]] auto cbegin() const -> EdgeConstIterator { return mEdges.cbegin(); }
  [[nodiscard]] auto cend() const -> EdgeConstIterator { return mEdges.cend(); }

 private:
  /// InlinedVector<u32, 2> stores up to 2 counts inline (no heap) —
  /// covers the standard 2-sample case. Spills to heap for >2 samples.
  using Counts = absl::InlinedVector<u32, 2>;

  // ── 8B Align ────────────────────────────────────────────────────────────
  EdgeList mEdges;
  Kmer mKmer;
  usize mCompId = 0;
  Counts mCounts;                    // 8B (InlinedVector<u32, 2>)
  std::array<u32, 2> mRoleCounts{};  // 8B (2×4B, [0]=CTRL [1]=CASE)

  // ── 1B Align ────────────────────────────────────────────────────────────
  Label mLabel;  // 1B
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_NODE_H_
