#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>

#include "absl/container/flat_hash_set.h"
#include "lancet/edge.h"
#include "lancet/graph.h"
#include "lancet/node.h"
#include "lancet/path.h"

namespace lancet {
class EdmondKarpMaxFlow {
 public:
  explicit EdmondKarpMaxFlow(const Graph::NodeContainer* nc, std::size_t kmer_size, std::size_t max_path_len,
                             std::uint32_t bfs_limit, bool is_tenx_mode = false);

  EdmondKarpMaxFlow() = delete;

  [[nodiscard]] auto NextPath() -> std::unique_ptr<Path>;

 private:
  const Graph::NodeContainer* nodesMap = nullptr;
  const Node* sourcePtr = nullptr;
  const Node* sinkPtr = nullptr;
  std::size_t kmerSize = 0;
  std::size_t maxPathLen = std::numeric_limits<std::size_t>::max();
  std::uint32_t bfsLimit = 0;
  bool isTenxMode = false;

  // to identify edges already returned in the previous call to `next_path`.
  absl::flat_hash_set<Edge> markedEdges;
};
}  // namespace lancet
