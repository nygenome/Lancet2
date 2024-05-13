#ifndef SRC_LANCET_CALLER_MSA_BUILDER_H_
#define SRC_LANCET_CALLER_MSA_BUILDER_H_

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "spoa/graph.hpp"

namespace lancet::caller {

class MsaBuilder {
 public:
  using FsPath = std::filesystem::path;
  using RefAndAltHaplotypes = absl::Span<const std::string>;

  explicit MsaBuilder(RefAndAltHaplotypes sequences, const FsPath& out_gfa_path = FsPath());

  [[nodiscard]] auto MultipleSequenceAlignment() const -> std::vector<std::string_view>;
  [[nodiscard]] auto FetchHaplotypeSeqView(const usize idx) const -> std::string_view { return mHaplotypeSeqs.at(idx); }

 private:
  RefAndAltHaplotypes mHaplotypeSeqs;
  std::vector<std::string> mResultMsa;

  static void WriteFasta(const FsPath& out_path, absl::Span<const std::string> msa_alns);
  static void WriteGfa(const FsPath& out_path, const spoa::Graph& graph);
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_MSA_BUILDER_H_
