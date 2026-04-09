#ifndef SRC_LANCET_CALLER_MSA_BUILDER_H_
#define SRC_LANCET_CALLER_MSA_BUILDER_H_

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"

namespace lancet::caller {

class MsaBuilder {
 public:
  std::unique_ptr<spoa::AlignmentEngine> mEngine;
  spoa::Graph mGraph;

  using FsPath = std::filesystem::path;

  void UpdateSpoaState(absl::Span<const std::string> sequences);
  void SerializeGraph(const FsPath& out_gfa_path);

 private:
  static void WriteFasta(const FsPath& out_path, absl::Span<const std::string> msa_alns);
  void WriteGfa(const FsPath& out_path) const;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_MSA_BUILDER_H_
