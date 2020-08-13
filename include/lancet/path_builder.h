#pragma once

#include <cstddef>
#include <memory>
#include <string>

#include "absl/container/inlined_vector.h"
#include "absl/types/span.h"
#include "lancet/core_enums.h"
#include "lancet/edge.h"
#include "lancet/node.h"
#include "lancet/path.h"

namespace lancet {
class PathBuilder {
 public:
  explicit PathBuilder(std::size_t k, bool is_tenx_mode);
  PathBuilder() = delete;

  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return nodesList.empty(); }
  [[nodiscard]] auto NumNodes() const noexcept -> std::size_t { return nodesList.size(); }
  [[nodiscard]] auto PathLength() const noexcept -> std::size_t { return pathLen; }
  [[nodiscard]] auto Direction() const noexcept -> Strand { return pathDir; }

  void MarkSinkTouch() { touchedSink = true; }
  [[nodiscard]] auto TouchedSink() const noexcept -> bool { return touchedSink; }

  [[nodiscard]] auto Score() const noexcept -> std::size_t { return pathScore; }
  void IncrementScore() { ++pathScore; }

  [[nodiscard]] auto LastNode() const noexcept -> const Node* {
    return nodesList.empty() ? nullptr : nodesList[nodesList.size() - 1];
  }

  void Extend(const Edge* link, const Node* destination);

  [[nodiscard]] auto PathEdges() -> absl::Span<const Edge*> { return absl::MakeSpan(edgesList); }
  [[nodiscard]] auto BuildPath() const -> std::unique_ptr<Path>;

 private:
  absl::InlinedVector<const Node*, 64> nodesList;
  absl::InlinedVector<const Edge*, 64> edgesList;
  Strand pathDir = Strand::FWD;
  std::size_t kmerSize = 0;
  std::size_t pathLen = 0;
  std::size_t pathScore = 0;
  bool isTenxMode = false;
  bool touchedSink = false;

  [[nodiscard]] auto BuildPathSeq() const -> std::string;
  [[nodiscard]] auto BuildPathCov() const -> NodeCov;
  [[nodiscard]] auto BuildPathHP() const -> NodeHP;
};
}  // namespace lancet
