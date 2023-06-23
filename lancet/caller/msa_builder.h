#ifndef SRC_LANCET_CALLER_MSA_BUILDER_H_
#define SRC_LANCET_CALLER_MSA_BUILDER_H_

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/hts/reference.h"
#include "spoa/spoa.hpp"

namespace lancet::caller {

class MsaBuilder {
 public:
  using FsPath = std::filesystem::path;
  using RefHaplotype = std::shared_ptr<const hts::Reference::Region>;
  using AltHaplotypes = absl::Span<const std::string>;

  explicit MsaBuilder(RefHaplotype ref_seq, AltHaplotypes alt_seqs, const FsPath& out_gfa_path = FsPath());

  [[nodiscard]] auto MultipleSequenceAlignment() const -> std::vector<std::string_view>;
  [[nodiscard]] auto FetchHaplotypeSeqView(const usize idx) const -> std::string_view {
    return idx == 0 ? mRefRegion->SeqView() : mAlternateSeqs.at(idx - 1);
  }

 private:
  RefHaplotype mRefRegion;
  AltHaplotypes mAlternateSeqs;
  std::vector<std::string> mResultMsa;

  static void WriteFasta(const FsPath& out_path, absl::Span<const std::string> msa_alns);
  static void WriteGfa(const FsPath& out_path, const spoa::Graph& graph, absl::Span<const bool> is_reversed);

  static void AddSequenceToGraph(const char* sequence, usize length, spoa::AlignmentEngine* aln_engine,
                                 spoa::Graph& graph, bool& added_reverse_complement);
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_MSA_BUILDER_H_
