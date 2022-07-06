#pragma once

#include <cstddef>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>

#include "absl/types/span.h"
#include "lancet2/graph.h"

namespace lancet2 {
class Graph::DotSerializer {
 public:
  explicit DotSerializer(const Graph* g) : graphPtr(g) {}

  void write_component(std::size_t comp_id, const std::string& suffix) const;
  void write_component(std::size_t comp_id, absl::Span<const PathNodeIds> flow_paths) const;

 private:
  const Graph* graphPtr = nullptr;

  static void dump_header(std::ostream& out_stream) ;
  void dump_component(std::size_t comp_id, std::ostream& out_stream) const;
  static void dump_path_flow(const PathNodeIds& path_flow, double hue, std::ostream& out_stream) ;

  static auto opposite_strand_sequence(std::string_view seq) -> std::string;
};
}  // namespace lancet2
