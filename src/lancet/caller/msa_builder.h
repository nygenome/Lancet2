#ifndef SRC_LANCET_CALLER_MSA_BUILDER_H_
#define SRC_LANCET_CALLER_MSA_BUILDER_H_

#include "lancet/base/types.h"

#include "absl/types/span.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace lancet::caller {

class MsaBuilder {
 public:
  std::unique_ptr<spoa::AlignmentEngine> mEngine;
  spoa::Graph mGraph;

  using FsPath = std::filesystem::path;

  void UpdateSpoaState(absl::Span<std::string const> sequences);
  void SerializeGraph(FsPath const& out_gfa_path);

 private:
  static void WriteFasta(FsPath const& gfa_path, absl::Span<std::string const> msa_alns);
  void WriteGfa(FsPath const& out_path) const;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_MSA_BUILDER_H_
