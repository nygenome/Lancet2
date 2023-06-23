#include "lancet/cbdg/max_flow.h"

#include <algorithm>

#include "absl/container/flat_hash_map.h"
#include "absl/hash/hash.h"
#include "absl/strings/cord.h"
#include "lancet/base/assert.h"

namespace lancet::cbdg {

MaxFlow::MaxFlow(const Graph::NodeTable *graph, const NodeIDPair &src_and_snk, const usize currk)
    : mGraph(graph), mCurrentK(currk) {
  const auto [source_id, sink_id] = src_and_snk;
  const auto src_itr = mGraph->find(source_id);
  LANCET_ASSERT(src_itr != mGraph->end())
  LANCET_ASSERT(src_itr->second != nullptr)
  mSource = src_itr->second.get();

  const auto snk_itr = mGraph->find(sink_id);
  LANCET_ASSERT(snk_itr != mGraph->end());
  LANCET_ASSERT(snk_itr->second != nullptr)
  mSink = snk_itr->second.get();
}

auto MaxFlow::NextPath() -> Result {
  auto walk = BuildNextWalk();
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (!walk.has_value() || walk->empty()) return std::nullopt;
  return BuildSequence(*walk);
}

auto MaxFlow::BuildNextWalk() -> std::optional<Walk> {
  using WalkID = usize;

  usize nvisits = 0;
  Walk best_possible_walk;
  CandidateWalks candidates;
  absl::flat_hash_map<WalkID, u64> walk_scores;

  static const auto make_walk_id = absl::Hash<Walk>();
  const auto dflt_src_sign = mSource->SignFor(Kmer::Ordering::DEFAULT);
  static constexpr usize ESTIMATED_WALK_LENGTH = 128;

  // Add outgoing edges from source node as seed walks for traversal
  PopulateWalkableEdgesInDirection(mSource, dflt_src_sign);
  for (const Edge &conn : mWalkableEdges) {
    Walk seed_walk;
    seed_walk.reserve(ESTIMATED_WALK_LENGTH);
    seed_walk.emplace_back(conn);
    candidates.emplace_back(std::move(seed_walk));

    const auto is_uniq_edge = mTraversed.find(conn) == mTraversed.cend();
    const auto identifier = make_walk_id(candidates.back());
    walk_scores.emplace(identifier, is_uniq_edge ? 1 : 0);
  }

  while (!candidates.empty()) {
    nvisits++;

    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (nvisits > Graph::DEFAULT_GRAPH_TRAVERSAL_LIMIT) break;

    const Walk &current_walk = candidates.front();
    const Edge &last_edge = current_walk.back();
    const Node *leaf_node = mGraph->at(last_edge.DstId()).get();

    const auto walk_direction = last_edge.DstSign();
    const auto candidate_walk_id = make_walk_id(current_walk);
    const auto current_score = walk_scores.at(candidate_walk_id);

    // If we touched sink node and have at-least one unique edge,
    // we set the current path as the best possible path and return
    if (leaf_node->Identifier() == mSink->Identifier() && current_score > 0) {
      walk_scores.erase(candidate_walk_id);
      best_possible_walk = candidates.front();
      break;
    }

    // If we touched sink node and do not have any unique edge,
    // we pop the current path from deque and keep searching
    if (leaf_node->Identifier() == mSink->Identifier() && current_score == 0) {
      walk_scores.erase(candidate_walk_id);
      candidates.pop_front();
      continue;
    }

    PopulateWalkableEdgesInDirection(leaf_node, walk_direction);
    for (const Edge &conn : mWalkableEdges) {
      Walk extension = current_walk;
      extension.reserve(ESTIMATED_WALK_LENGTH);
      extension.emplace_back(conn);
      candidates.emplace_back(std::move(extension));

      const auto is_uniq_edge = mTraversed.find(conn) == mTraversed.end();
      const auto identifier = make_walk_id(candidates.back());
      walk_scores.emplace(identifier, is_uniq_edge ? current_score + 1 : current_score);
    }

    LANCET_ASSERT(!candidates.empty())
    walk_scores.erase(candidate_walk_id);
    candidates.pop_front();
  }

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (best_possible_walk.empty()) return std::nullopt;

  mTraversed.insert(best_possible_walk.cbegin(), best_possible_walk.cend());
  return best_possible_walk;
}

auto MaxFlow::BuildSequence(WalkView walk) const -> Result {
  absl::Cord merged_seq;
  LANCET_ASSERT(!walk.empty())

  constexpr auto dflt_order = Kmer::Ordering::DEFAULT;
  constexpr auto oppo_order = Kmer::Ordering::OPPOSITE;
  auto ordering = walk[0].SrcSign() == Kmer::Sign::PLUS ? dflt_order : oppo_order;

  for (const auto &conn : walk) {
    if (merged_seq.empty()) {
      const auto src_itr = mGraph->find(conn.SrcId());
      LANCET_ASSERT(src_itr != mGraph->end())
      LANCET_ASSERT(src_itr->second != nullptr)
      merged_seq.Append(src_itr->second->CordDataFor(ordering));
    }

    const auto dst_itr = mGraph->find(conn.DstId());
    LANCET_ASSERT(dst_itr != mGraph->end())
    LANCET_ASSERT(dst_itr->second != nullptr)

    ordering = conn.DstSign() == Kmer::Sign::PLUS ? dflt_order : oppo_order;
    const auto dst_cord = dst_itr->second->CordDataFor(ordering);
    const auto uniq_seq_len = dst_cord.size() - mCurrentK + 1;
    merged_seq.Append(dst_cord.Subcord(mCurrentK - 1, uniq_seq_len));
  }

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (merged_seq.empty()) return std::nullopt;

  return std::string(merged_seq);
}

void MaxFlow::PopulateWalkableEdgesInDirection(const Node *src, Kmer::Sign dir) {
  LANCET_ASSERT(src != nullptr)
  mWalkableEdges.clear();
  mWalkableEdges.reserve(src->NumOutEdges());

  std::ranges::copy_if(*src, std::back_inserter(mWalkableEdges),
                       [&dir](const Edge &edge) -> bool { return edge.SrcSign() == dir; });

  std::ranges::sort(mWalkableEdges, [this](const Edge &lhs, const Edge &rhs) -> bool {
    // Extend candidate paths with bfs. Sort node edges by prioritizing
    // unwalked paths first and then sort later for deterministic traversal
    const auto is_unwalked_lhs = this->mTraversed.find(lhs) == this->mTraversed.end();
    const auto is_unwalked_rhs = this->mTraversed.find(rhs) == this->mTraversed.end();
    // NOLINTBEGIN(readability-braces-around-statements)
    if (is_unwalked_lhs != is_unwalked_rhs) return is_unwalked_lhs;
    if (lhs.SrcId() != rhs.SrcId()) return lhs.SrcId() < rhs.SrcId();
    if (lhs.DstId() != rhs.DstId()) return lhs.DstId() < rhs.DstId();
    return static_cast<u64>(lhs.Kind()) < static_cast<u64>(rhs.Kind());
    // NOLINTEND(readability-braces-around-statements)
  });
}

}  // namespace lancet::cbdg
