#pragma once

#include <memory>
#include <string>

#include "absl/container/inlined_vector.h"
#include "absl/types/span.h"
#include "lancet2/core_enums.h"
#include "lancet2/edge.h"
#include "lancet2/node.h"
#include "lancet2/path.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class PathBuilder {
 public:
  explicit PathBuilder(usize k, bool is_tenx_mode);
  PathBuilder() = delete;

  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return nodesList.empty(); }
  [[nodiscard]] auto NumNodes() const noexcept -> usize { return nodesList.size(); }
  [[nodiscard]] auto PathLength() const noexcept -> usize { return pathLen; }
  [[nodiscard]] auto Direction() const noexcept -> Strand { return pathDir; }

  void MarkSinkTouch() { touchedSink = true; }
  [[nodiscard]] auto TouchedSink() const noexcept -> bool { return touchedSink; }

  [[nodiscard]] auto Score() const noexcept -> usize { return pathScore; }
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
  usize kmerSize = 0;
  usize pathLen = 0;
  usize pathScore = 0;
  bool isTenxMode = false;
  bool touchedSink = false;

  [[nodiscard]] auto BuildPathSeq() const -> std::string;
};
}  // namespace lancet2
