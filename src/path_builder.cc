#include "lancet2/path_builder.h"

#include <utility>

#include "lancet2/assert_macro.h"
#include "lancet2/merge_node_info.h"
#include "lancet2/utils.h"

namespace lancet2 {
PathBuilder::PathBuilder(std::size_t k, bool is_tenx_mode) : kmerSize(k), isTenxMode(is_tenx_mode) {}

void PathBuilder::Extend(const Edge *link, const Node *destination) {
  LANCET_ASSERT(!destination->IsMockNode());  // NOLINT
  nodesList.push_back(destination);
  edgesList.push_back(link);
  pathDir = link->DstDirection();
  pathLen += pathLen == 0 ? destination->Length() : (destination->Length() - kmerSize + 1);
}

auto PathBuilder::BuildPath() const -> std::unique_ptr<Path> {
  auto pathSeq = BuildPathSeq();
  if (pathSeq.empty()) return nullptr;
  return std::make_unique<Path>(absl::FixedArray<const Node *>(nodesList.cbegin(), nodesList.cend()),
                                absl::FixedArray<const Edge *>(edgesList.cbegin(), edgesList.cend()),
                                std::move(pathSeq), BuildPathCov(), isTenxMode ? BuildPathHP() : NodeHP{});
}

auto PathBuilder::BuildPathSeq() const -> std::string {
  std::string result;
  result.reserve(pathLen);
  LANCET_ASSERT(!nodesList.empty() && !edgesList.empty() && nodesList.size() == edgesList.size());  // NOLINT

  for (const auto &node : nodesList) {
    if (result.empty()) {
      node->Orientation() == Strand::REV ? utils::PushRevCompSeq(node->SeqView(), &result)
                                         : utils::PushSeq(node->SeqView(), &result);
      continue;
    }

    const auto shouldRev = node->Orientation() == Strand::REV;
    if (!CanMergeSeqs(result, node->SeqView(), BuddyPosition::FRONT, shouldRev, kmerSize)) {
      // Found irrecoverable error/cycle in path. break and return empty path
      return {};
    }

    shouldRev ? utils::PushRevCompSeq(&result, result.end(), node->SeqView(), kmerSize - 1, node->Length())
              : utils::PushSeq(&result, result.end(), node->SeqView(), kmerSize - 1, node->Length());
  }

  LANCET_ASSERT(result.length() == pathLen);  // NOLINT
  return result;
}

auto PathBuilder::BuildPathCov() const -> NodeCov {
  NodeCov result;
  result.Reserve(pathLen);
  for (const auto &node : nodesList) {
    result.MergeBuddy(node->CovData(), BuddyPosition::FRONT, node->Orientation() == Strand::REV, kmerSize);
  }
  LANCET_ASSERT(result.Size() == pathLen);  // NOLINT
  return result;
}

auto PathBuilder::BuildPathHP() const -> NodeHP {
  NodeHP result;
  result.Reserve(pathLen);
  for (const auto &node : nodesList) {
    result.MergeBuddy(node->HPData(), BuddyPosition::FRONT, node->Orientation() == Strand::REV, kmerSize);
  }
  LANCET_ASSERT(result.Size() == pathLen);  // NOLINT
  return result;
}
}  // namespace lancet2
