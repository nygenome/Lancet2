#include "lancet2/path_builder.h"

#include <utility>

#include "lancet2/assert_macro.h"
#include "lancet2/merge_node_info.h"
#include "lancet2/utils.h"

namespace lancet2 {
PathBuilder::PathBuilder(usize k, bool is_tenx_mode) : kmerSize(k), isTenxMode(is_tenx_mode) {}

void PathBuilder::Extend(const Edge *link, const Node *destination) {
  LANCET_ASSERT(!destination->IsMockNode());  // NOLINT
  nodesList.push_back(destination);
  edgesList.push_back(link);
  pathDir = link->GetDstDir();
  pathLen += pathLen == 0 ? destination->GetLength() : (destination->GetLength() - kmerSize + 1);
}

auto PathBuilder::BuildPath() const -> std::unique_ptr<Path> {
  auto pathSeq = BuildPathSeq();
  if (pathSeq.empty()) return nullptr;
  return std::make_unique<Path>(absl::FixedArray<const Node *>(nodesList.cbegin(), nodesList.cend()),
                                absl::FixedArray<const Edge *>(edgesList.cbegin(), edgesList.cend()),
                                std::move(pathSeq));
}

auto PathBuilder::BuildPathSeq() const -> std::string {
  std::string result;
  result.reserve(pathLen);
  LANCET_ASSERT(!nodesList.empty() && !edgesList.empty() && nodesList.size() == edgesList.size());  // NOLINT

  for (const auto &node : nodesList) {
    if (result.empty()) {
      node->GetOrientation() == Strand::REV ? utils::PushRevCompSeq(node->GetSeqView(), &result)
                                            : utils::PushSeq(node->GetSeqView(), &result);
      continue;
    }

    const auto shouldRev = node->GetOrientation() == Strand::REV;
    if (!CanMergeSeqs(result, node->GetSeqView(), BuddyPosition::FRONT, shouldRev, kmerSize)) {
      // Found irrecoverable error/cycle in path. break and return empty path
      return {};
    }

    shouldRev ? utils::PushRevCompSeq(&result, result.end(), node->GetSeqView(), kmerSize - 1, node->GetLength())
              : utils::PushSeq(&result, result.end(), node->GetSeqView(), kmerSize - 1, node->GetLength());
  }

  LANCET_ASSERT(result.length() == pathLen);  // NOLINT
  return result;
}
}  // namespace lancet2
