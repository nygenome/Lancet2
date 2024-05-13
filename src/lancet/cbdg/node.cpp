#include "lancet/cbdg/node.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <vector>

#include "lancet/base/types.h"
#include "lancet/cbdg/kmer.h"

namespace {

inline auto WeightedAverage(const std::array<u32, 2>& numbers, const std::array<usize, 2>& weights) -> u32 {
  return ((numbers[0] * weights[0]) + (numbers[1] * weights[1])) / (weights[0] + weights[1]);
}

}  // namespace

namespace lancet::cbdg {

void Node::AddLabel(const Label& label) { mLabel.Merge(label); }

void Node::IncrementReadSupport(const Label& label) {
  if (label.HasTag(Label::NORMAL)) {
    mCounts[NORMAL_COUNT_INDEX] += 1;
  }

  if (label.HasTag(Label::TUMOR)) {
    mCounts[TUMOR_COUNT_INDEX] += 1;
  }
}

auto Node::NormalReadSupport() const noexcept -> u32 { return mCounts[NORMAL_COUNT_INDEX]; }
auto Node::TumorReadSupport() const noexcept -> u32 { return mCounts[TUMOR_COUNT_INDEX]; }
auto Node::TotalReadSupport() const noexcept -> u32 { return NormalReadSupport() + TumorReadSupport(); }

void Node::Merge(const Node& other, const EdgeKind conn_kind, const usize currk) {
  mKmer.Merge(other.mKmer, conn_kind, currk);
  mLabel.Merge(other.mLabel);
  mCounts[0] = WeightedAverage({mCounts[0], other.mCounts[0]}, {mKmer.Length(), other.Length()});
  mCounts[1] = WeightedAverage({mCounts[1], other.mCounts[1]}, {mKmer.Length(), other.Length()});
}

auto Node::HasSelfLoop() const -> bool {
  return std::ranges::any_of(mEdges, [](const Edge& conn) -> bool { return conn.IsSelfLoop(); });
}

auto Node::FindEdgesInDirection(const Kmer::Ordering ord) const -> std::vector<Edge> {
  std::vector<Edge> results;
  results.reserve(mEdges.size());
  const auto expected_src_sign = mKmer.SignFor(ord);
  std::ranges::copy_if(mEdges, std::back_inserter(results),
                       [&expected_src_sign](const Edge& conn) -> bool { return conn.SrcSign() == expected_src_sign; });
  return results;
}

}  // namespace lancet::cbdg
