#pragma once

#include <memory>

#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "lancet2/base_hpcov.h"
#include "lancet2/cli_params.h"
#include "lancet2/graph.h"
#include "lancet2/node.h"
#include "lancet2/read_info.h"
#include "lancet2/ref_window.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class GraphBuilder {
 public:
  GraphBuilder(std::shared_ptr<const RefWindow> w, absl::Span<const ReadInfo> reads, double avg_cov,
               std::shared_ptr<const CliParams> p);
  GraphBuilder() = delete;

  [[nodiscard]] auto BuildGraph(usize min_k, usize max_k) -> std::unique_ptr<Graph>;
  [[nodiscard]] auto CurrentKmerSize() const noexcept -> usize { return currentK; }

 private:
  double avgCov = 0.0;
  usize currentK = 0;
  std::shared_ptr<const RefWindow> window;
  std::shared_ptr<const CliParams> params;
  absl::Span<const ReadInfo> sampleReads;
  Graph::NodeContainer nodesMap;

  void BuildSampleNodes();
  void BuildRefNodes();

  void RecoverKmers();

  [[nodiscard]] auto BuildNodes(absl::string_view seq) -> std::vector<NodeIdentifier>;

  void BuildNode(NodeIdentifier node_id);

  [[nodiscard]] static auto MutateSeq(absl::string_view seq, usize base_pos) -> std::vector<std::string>;
};
}  // namespace lancet2
