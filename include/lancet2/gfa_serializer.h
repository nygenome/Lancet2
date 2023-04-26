#pragma once

#include <ostream>
#include <string>
#include <string_view>
#include <utility>

#include "absl/types/span.h"
#include "lancet2/graph.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class Graph::GfaSerializer {
 public:
  explicit GfaSerializer(const Graph* g) : graphPtr(g) {}

  void WriteComponent(usize comp_id, const std::string& suffix) const;
  void WriteComponent(usize comp_id, absl::Span<const PathNodeIds> flow_paths, std::string& windowId) const;

 private:
  const Graph* graphPtr = nullptr;

  static void DumpHeader(std::ostream& out_stream);
  void DumpComponent(usize comp_id, std::ostream& out_stream) const;
  void DumpPathFlow(const PathNodeIds& path_flow, usize pathNum, std::string& windowId, std::ostream& out_stream) const;

  static auto OppositeStrandSequence(std::string_view seq) -> std::string;
};
}  // namespace lancet2
