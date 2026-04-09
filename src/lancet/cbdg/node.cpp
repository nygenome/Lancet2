#include "lancet/cbdg/node.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/label.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <vector>

namespace {

inline auto WeightedAverage(std::array<u32, 2> const& numbers, std::array<usize, 2> const& weights)
    -> u32 {
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
  return ((numbers[0] * weights[0]) + (numbers[1] * weights[1])) / (weights[0] + weights[1]);
}

}  // namespace

namespace lancet::cbdg {

void Node::AddLabel(Label const& label) {
  mLabel.Merge(label);
}

void Node::IncrementReadSupport(Label const& label) {
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
  if (label.HasTag(Label::NORMAL)) {
    mCounts[NORMAL_COUNT_INDEX] += 1;
  }

  if (label.HasTag(Label::TUMOR)) {
    mCounts[TUMOR_COUNT_INDEX] += 1;
  }
}

auto Node::NormalReadSupport() const noexcept -> u32 {
  return mCounts[NORMAL_COUNT_INDEX];
}
auto Node::TumorReadSupport() const noexcept -> u32 {
  return mCounts[TUMOR_COUNT_INDEX];
}
// NOLINTEND(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
auto Node::TotalReadSupport() const noexcept -> u32 {
  return NormalReadSupport() + TumorReadSupport();
}

void Node::Merge(Node const& other, EdgeKind const conn_kind, usize const currk) {
  mKmer.Merge(other.mKmer, conn_kind, currk);
  mLabel.Merge(other.mLabel);
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
  mCounts[0] = WeightedAverage({mCounts[0], other.mCounts[0]}, {mKmer.Length(), other.Length()});
  mCounts[1] = WeightedAverage({mCounts[1], other.mCounts[1]}, {mKmer.Length(), other.Length()});
  // NOLINTEND(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
}

auto Node::HasSelfLoop() const -> bool {
  return std::ranges::any_of(mEdges, [](Edge const& conn) -> bool { return conn.IsSelfLoop(); });
}

auto Node::FindEdgesInDirection(Kmer::Ordering const ord) const -> std::vector<Edge> {
  std::vector<Edge> results;
  results.reserve(mEdges.size());
  auto const expected_src_sign = mKmer.SignFor(ord);
  std::ranges::copy_if(mEdges, std::back_inserter(results),
                       [&expected_src_sign](Edge const& conn) -> bool {
                         return conn.SrcSign() == expected_src_sign;
                       });
  return results;
}

}  // namespace lancet::cbdg
