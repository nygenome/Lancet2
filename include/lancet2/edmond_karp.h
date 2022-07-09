#pragma once

#include <limits>

#include "absl/container/flat_hash_set.h"
#include "lancet2/edge.h"
#include "lancet2/graph.h"
#include "lancet2/node.h"
#include "lancet2/path.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class EdmondKarpMaxFlow {
 public:
  explicit EdmondKarpMaxFlow(const Graph::NodeContainer* nc, usize kmer_size, usize max_path_len, u32 bfs_limit,
                             bool is_tenx_mode = false);

  EdmondKarpMaxFlow() = delete;

  [[nodiscard]] auto NextPath() -> std::unique_ptr<Path>;

 private:
  const Graph::NodeContainer* nodesMap = nullptr;
  const Node* sourcePtr = nullptr;
  usize kmerSize = 0;
  usize maxPathLen = std::numeric_limits<usize>::max();
  u32 bfsLimit = 0;
  bool isTenxMode = false;

  // to identify edges already returned to the previous call to `next_path`.
  absl::flat_hash_set<const Edge*> markedEdges;
};
}  // namespace lancet2
