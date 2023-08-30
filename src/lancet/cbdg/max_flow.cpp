#include "lancet/cbdg/max_flow.h"

#include <algorithm>

#include "absl/container/flat_hash_map.h"
#include "absl/hash/hash.h"
#include "absl/strings/str_cat.h"
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
  LANCET_ASSERT(!walk.empty())

  usize total_seq_len = 0;
  std::vector<std::string> uniq_seqs;
  uniq_seqs.reserve(walk.size() + 1);

  constexpr auto dflt_order = Kmer::Ordering::DEFAULT;
  constexpr auto oppo_order = Kmer::Ordering::OPPOSITE;
  auto ordering = walk[0].SrcSign() == Kmer::Sign::PLUS ? dflt_order : oppo_order;

  for (const auto &conn : walk) {
    if (uniq_seqs.empty()) {
      const auto src_itr = mGraph->find(conn.SrcId());
      LANCET_ASSERT(src_itr != mGraph->end())
      LANCET_ASSERT(src_itr->second != nullptr)
      uniq_seqs.emplace_back(src_itr->second->SequenceFor(ordering));
      total_seq_len += uniq_seqs.back().length();
    }

    const auto dst_itr = mGraph->find(conn.DstId());
    LANCET_ASSERT(dst_itr != mGraph->end())
    LANCET_ASSERT(dst_itr->second != nullptr)

    ordering = conn.DstSign() == Kmer::Sign::PLUS ? dflt_order : oppo_order;
    const auto dst_seq = dst_itr->second->SequenceFor(ordering);
    const auto uniq_seq_len = dst_seq.size() - mCurrentK + 1;
    uniq_seqs.emplace_back(dst_seq.substr(mCurrentK - 1, uniq_seq_len));
    total_seq_len += uniq_seqs.back().length();
  }

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (uniq_seqs.empty()) return std::nullopt;

  std::string merged_seq;
  merged_seq.reserve(total_seq_len);
  for (const auto &item : uniq_seqs) {
    absl::StrAppend(&merged_seq, item);
  }

  LANCET_ASSERT(merged_seq.length() == total_seq_len)
  return merged_seq;
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
