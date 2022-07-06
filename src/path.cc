#include "lancet2/path.h"

#include <algorithm>
#include <utility>

#include "lancet2/assert_macro.h"

namespace lancet2 {
Path::Path(absl::FixedArray<const Node *> path_nodes, absl::FixedArray<const Edge *> path_edges, std::string path_seq,
           NodeCov path_covs, NodeHP path_hps)
    : nodesList(std::move(path_nodes)), edgesList(std::move(path_edges)), pathSeq(std::move(path_seq)),
      pathCovs(std::move(path_covs)), pathHPs(std::move(path_hps)) {}

auto Path::FindSpanningNode(usize path_pos, usize curr_k) const -> const Node * {
  const auto *result = std::find_if(nodesList.cbegin(), nodesList.cend(), [&path_pos, &curr_k](const auto &node) {
    static usize currPos = 0;
    if (node->IsMockNode()) return false;
    if (currPos + node->Length() >= path_pos) return true;
    currPos += node->Length() - curr_k + 1;
    return false;
  });

  return result == nodesList.cend() ? nullptr : *result;
}

auto Path::TouchedEdgeIDs() const -> PathNodeIds {
  absl::FixedArray<EdgeNodeIds> result(edgesList.size() + 1);

  LANCET_ASSERT(!nodesList.empty() && !edgesList.empty() && nodesList.size() == edgesList.size());  // NOLINT

  auto currentSrcId = MOCK_SOURCE_ID;
  for (usize idx = 0; idx < nodesList.size(); idx++) {
    result[idx].srcId = currentSrcId;
    result[idx].dstId = edgesList[idx]->DestinationID();
    currentSrcId = edgesList[idx]->DestinationID();
  }

  result[result.size() - 1].srcId = currentSrcId;
  result[result.size() - 1].dstId = MOCK_SINK_ID;

  return result;
}
}  // namespace lancet2
