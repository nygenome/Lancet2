#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "lancet/base_hpcov.h"
#include "lancet/cli_params.h"
#include "lancet/graph.h"
#include "lancet/node.h"
#include "lancet/read_info.h"
#include "lancet/ref_window.h"

namespace lancet {
class GraphBuilder {
 public:
  GraphBuilder(std::shared_ptr<const RefWindow> w, absl::Span<const ReadInfo> reads, double avg_cov,
               std::shared_ptr<const CliParams> p);
  GraphBuilder() = delete;

  [[nodiscard]] auto BuildGraph(std::size_t min_k, std::size_t max_k) -> std::unique_ptr<Graph>;

  using ReferenceData = std::vector<BaseHpCov>;
  [[nodiscard]] auto RefData(SampleLabel label) const noexcept -> ReferenceData {
    return label == SampleLabel::NORMAL ? refNmlData : refTmrData;
  }

  [[nodiscard]] auto CurrentKmerSize() const noexcept -> std::size_t { return currentK; }

 private:
  std::shared_ptr<const RefWindow> window;
  absl::Span<const ReadInfo> sampleReads;
  std::shared_ptr<const CliParams> params;
  double avgCov = 0.0;
  std::size_t currentK = 0;
  Graph::NodeContainer nodesMap;
  ReferenceData refTmrData;
  ReferenceData refNmlData;

  void BuildSampleNodes();
  void BuildRefNodes();

  void RecoverKmers();

  struct BuildNodesResult {
    std::vector<NodeIdentifier> nodeIDs;
    std::size_t numNodesBuilt = 0;
    std::size_t numKmersGiven = 0;
  };
  [[nodiscard]] auto BuildNodes(absl::string_view seq) -> BuildNodesResult;

  struct BuildNodeResult {
    NodeIdentifier ID = 0;  // NOLINT
    bool builtNode = false;
  };
  auto BuildNode(NodeIdentifier node_id) -> BuildNodeResult;

  void BuildRefData(absl::Span<const std::size_t> ref_mer_hashes);

  [[nodiscard]] static auto MutateSeq(absl::string_view seq, std::size_t base_pos) -> std::vector<std::string>;
};
}  // namespace lancet
