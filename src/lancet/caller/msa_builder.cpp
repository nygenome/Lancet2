#include "lancet/caller/msa_builder.h"

#include <algorithm>
#include <fstream>
#include <utility>

#include "spdlog/fmt/fmt.h"
#include "spdlog/fmt/ostr.h"

namespace {

using AlnEngine = std::unique_ptr<spoa::AlignmentEngine>;
inline auto BuildAlnEngine(absl::Span<const std::string> seqs) -> AlnEngine {
  static constexpr i8 MATCH = 1;
  static constexpr i8 MISMATCH = -19;
  static constexpr i8 OPEN1 = -81;
  static constexpr i8 EXTEND1 = -1;
  static constexpr i8 OPEN2 = -39;
  static constexpr i8 EXTEND2 = -3;
  // asm5 from minimap2 -> https://lh3.github.io/minimap2/minimap2.html -> assembly to same species ref scoring
  // https://curiouscoding.nl/posts/pairwise-alignment -> Convex affine gap scoring -> min(g1+(i-1)*e1, g2+(i-1)*e2)
  auto aln = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, MATCH, MISMATCH, OPEN1, EXTEND1, OPEN2, EXTEND2);
  const auto* longest_seq = std::ranges::max_element(seqs, std::less<>(), &std::string::size);
  aln->Prealloc(longest_seq->length(), 4);
  return aln;
}

}  // namespace

namespace lancet::caller {

MsaBuilder::MsaBuilder(RefAndAltHaplotypes sequences, const FsPath& out_gfa_path) : mHaplotypeSeqs(sequences) {
  mResultMsa.reserve(mHaplotypeSeqs.size());

  spoa::Graph graph;
  const auto engine = BuildAlnEngine(mHaplotypeSeqs);
  std::ranges::for_each(sequences, [&engine, &graph](const std::string& haplotype) {
    const auto alignment = engine->Align(haplotype, graph);
    graph.AddAlignment(alignment, haplotype);
  });

  mResultMsa = graph.GenerateMultipleSequenceAlignment(false);
  if (!out_gfa_path.empty()) {
    WriteFasta(out_gfa_path, mResultMsa);
    WriteGfa(out_gfa_path, graph);
  }
}

auto MsaBuilder::MultipleSequenceAlignment() const -> std::vector<std::string_view> {
  std::vector<std::string_view> results_view;
  results_view.reserve(mResultMsa.size());
  std::ranges::transform(mResultMsa, std::back_inserter(results_view),
                         [](const std::string& path) -> std::string_view { return path; });
  return results_view;
}

void MsaBuilder::WriteFasta(const FsPath& gfa_path, absl::Span<const std::string> msa_alns) {
  auto fa_path = gfa_path.parent_path() / gfa_path.stem();
  fa_path += ".fasta";
  std::ofstream out_handle(fa_path);
  for (usize idx = 0; idx < msa_alns.size(); ++idx) {
    fmt::print(out_handle, ">{}{}\n{}\n", idx == 0 ? "ref" : "hap", idx, msa_alns[idx]);
  }
  out_handle.close();
}

void MsaBuilder::WriteGfa(const FsPath& out_path, const spoa::Graph& graph) {
  // https://github.com/rvaser/spoa/pull/36
  // See PR for how to normalize & process the output GFA
  std::ofstream out_handle(out_path, std::ios::trunc);
  fmt::print(out_handle, "H\tVN:Z:1.0\n");

  for (const std::unique_ptr<spoa::Graph::Node>& node : graph.nodes()) {
    fmt::print(out_handle, "S\t{}\t{}\n", node->id + 1, static_cast<char>(graph.decoder(node->code)));
    for (const spoa::Graph::Edge* edge : node->outedges) {
      fmt::print(out_handle, "L\t{}\t+\t{}\t+\t0M\n", node->id + 1, edge->head->id + 1);
    }
  }

  for (u32 seq_idx = 0; seq_idx < graph.sequences().size(); ++seq_idx) {
    fmt::print(out_handle, "P\t{}{}\t", seq_idx == 0 ? "ref" : "hap", seq_idx);

    std::vector<std::uint32_t> path;
    const spoa::Graph::Node* current_node = graph.sequences()[seq_idx];
    while (current_node != nullptr) {
      path.emplace_back(current_node->id + 1);
      current_node = current_node->Successor(seq_idx);
    }

    for (usize path_idx = 0; path_idx < path.size(); ++path_idx) {
      fmt::print(out_handle, "{}{}+", path_idx == 0 ? "" : ",", path[path_idx]);
    }

    fmt::print(out_handle, "\t*\n");
  }

  out_handle.close();
}

}  // namespace lancet::caller
