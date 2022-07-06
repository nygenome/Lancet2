#include "lancet2/edmond_karp.h"

#include <deque>
#include <utility>

#include "lancet2/assert_macro.h"
#include "lancet2/core_enums.h"
#include "lancet2/log_macros.h"
#include "lancet2/path_builder.h"
#include "spdlog/spdlog.h"

namespace lancet2 {
EdmondKarpMaxFlow::EdmondKarpMaxFlow(const Graph::NodeContainer *nc, usize kmer_size, usize max_path_len, u32 bfs_limit,
                                     bool is_tenx_mode)
    : nodesMap(nc), kmerSize(kmer_size), maxPathLen(max_path_len), bfsLimit(bfs_limit), isTenxMode(is_tenx_mode) {
  LANCET_ASSERT(nodesMap != nullptr);  // NOLINT

  const auto srcItr = nodesMap->find(MOCK_SOURCE_ID);
  LANCET_ASSERT(srcItr != nodesMap->end() && srcItr->second != nullptr);  // NOLINT
  LANCET_ASSERT(srcItr->second->NumEdges() == 1);                         // NOLINT
  LANCET_ASSERT(srcItr->second->NumEdges(Strand::FWD) == 1);              // NOLINT
  sourcePtr = srcItr->second.get();

#ifndef NDEBUG
  const auto snkItr = nodesMap->find(MOCK_SINK_ID);
  LANCET_ASSERT(snkItr != nodesMap->end() && snkItr->second != nullptr);  // NOLINT
  LANCET_ASSERT(snkItr->second->NumEdges() == 1);                         // NOLINT
#endif
}

auto EdmondKarpMaxFlow::NextPath() -> std::unique_ptr<Path> {
  u32 numVisits = 0;
  PathBuilder bestBuilder(kmerSize, isTenxMode);
  std::deque<PathBuilder> candidateBuilders;
  candidateBuilders.emplace_back(kmerSize, isTenxMode);

  while (!candidateBuilders.empty()) {
    numVisits++;
    if (numVisits > bfsLimit) break;

    auto &currBuilder = candidateBuilders.front();
    const auto *lastNode = (currBuilder.NumNodes() == 0 && numVisits == 1) ? sourcePtr : currBuilder.LastNode();
    LANCET_ASSERT(lastNode != nullptr);  // NOLINT

    if (currBuilder.PathLength() > maxPathLen) {
      // we extended the path too long. we don't care anymore
      candidateBuilders.pop_front();
      continue;
    }

    if (currBuilder.TouchedSink() && currBuilder.Score() > 0) {
      bestBuilder = currBuilder;
      break;
    }

    for (const Edge &e : *lastNode) {
      if (e.DestinationID() == MOCK_SINK_ID) {
        if (currBuilder.Score() <= bestBuilder.Score()) continue;
        PathBuilder srcToSink(currBuilder);
        srcToSink.MarkSinkTouch();
        candidateBuilders.emplace_back(std::move(srcToSink));
        continue;
      }

      if (e.DestinationID() == MOCK_SOURCE_ID || e.SrcDirection() != currBuilder.Direction()) continue;
      const auto neighbourItr = nodesMap->find(e.DestinationID());
      LANCET_ASSERT(neighbourItr != nodesMap->end());

      PathBuilder extensionBuilder(currBuilder);
      // If extension found a new edge between multiple calls to next_path, increment path score
      const auto uniqEdgeTouched = markedEdges.find(&e) == markedEdges.end();
      if (uniqEdgeTouched) extensionBuilder.IncrementScore();

      extensionBuilder.Extend(&e, neighbourItr->second.get());
      candidateBuilders.emplace_back(std::move(extensionBuilder));
    }

    LANCET_ASSERT(!candidateBuilders.empty());
    candidateBuilders.pop_front();
  }

  LOG_TRACE("Exiting Edmond Karp traversal after {} visits", numVisits);
  if (bestBuilder.IsEmpty()) return nullptr;
  const auto bestPathEdges = bestBuilder.PathEdges();
  markedEdges.insert(bestPathEdges.cbegin(), bestPathEdges.cend());
  return bestBuilder.BuildPath();
}
}  // namespace lancet2
