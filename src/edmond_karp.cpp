#include "lancet/edmond_karp.h"

#include <cassert>
#include <deque>
#include <utility>

#include "lancet/core_enums.h"
#include "lancet/path_builder.h"

namespace lancet {
EdmondKarpMaxFlow::EdmondKarpMaxFlow(const Graph::NodeContainer *nc, std::size_t kmer_size, std::size_t max_path_len,
                                     std::uint32_t bfs_limit, bool is_tenx_mode)
    : nodesMap(nc), kmerSize(kmer_size), maxPathLen(max_path_len), bfsLimit(bfs_limit), isTenxMode(is_tenx_mode) {
  assert(nodesMap != nullptr);  // NOLINT

  const auto srcItr = nodesMap->find(MOCK_SOURCE_ID);
  assert(srcItr != nodesMap->end() && srcItr->second != nullptr);  // NOLINT
  assert(srcItr->second->NumEdges() == 1);                         // NOLINT
  assert(srcItr->second->NumEdges(Strand::FWD) == 1);              // NOLINT

  const auto snkItr = nodesMap->find(MOCK_SINK_ID);
  assert(snkItr != nodesMap->end() && snkItr->second != nullptr);  // NOLINT
  assert(snkItr->second->NumEdges() == 1);                         // NOLINT

  sourcePtr = srcItr->second.get();
  sinkPtr = snkItr->second.get();
}

auto EdmondKarpMaxFlow::NextPath() -> std::unique_ptr<Path> {
  std::uint32_t numVisits = 0;
  PathBuilder bestBuilder(kmerSize, isTenxMode);
  std::deque<PathBuilder> candidateBuilders;
  candidateBuilders.emplace_back(kmerSize, isTenxMode);

  while (!candidateBuilders.empty()) {
    numVisits++;
    if (numVisits > bfsLimit) break;

    auto &currBuilder = candidateBuilders.front();
    const auto *lastNode = (currBuilder.NumNodes() == 0 && numVisits == 1) ? sourcePtr : currBuilder.LastNode();
    assert(lastNode != nullptr);  // NOLINT

    if (currBuilder.TouchedSink() && currBuilder.Score() > 0) {
      if (bestBuilder.IsEmpty() || currBuilder.Score() > bestBuilder.Score()) bestBuilder = currBuilder;
    } else if (currBuilder.PathLength() > maxPathLen) {
      // we extended the path too long. we don't care anymore
      candidateBuilders.pop_front();
      continue;
    } else {
      for (const Edge &e : *lastNode) {
        if (e.DestinationID() == MOCK_SINK_ID) {
          PathBuilder srcToSink(currBuilder);
          srcToSink.MarkSinkTouch();
          candidateBuilders.emplace_back(std::move(srcToSink));
          continue;
        }

        if (e.DestinationID() == MOCK_SOURCE_ID || e.SrcDirection() != currBuilder.Direction()) continue;
        PathBuilder extensionBuilder(currBuilder);
        // If extension found a new edge between multiple calls to next_path, increment path score
        const auto uniqEdgeTouched = markedEdges.find(e) == markedEdges.end();
        if (uniqEdgeTouched) extensionBuilder.IncrementScore();

        const auto neighbourItr = nodesMap->find(e.DestinationID());
        assert(neighbourItr != nodesMap->end());  // NOLINT
        extensionBuilder.Extend(e, neighbourItr->second.get());
        candidateBuilders.emplace_back(std::move(extensionBuilder));
      }
    }

    candidateBuilders.pop_front();
  }

  if (bestBuilder.IsEmpty()) return nullptr;
  const auto bestPathEdges = bestBuilder.PathEdges();
  markedEdges.insert(bestPathEdges.cbegin(), bestPathEdges.cend());
  return bestBuilder.BuildPath();
}
}  // namespace lancet
