#include "lancet/caller/msa_builder.h"

#include <algorithm>
#include <fstream>
#include <utility>

#include "absl/container/fixed_array.h"
#include "lancet/base/rev_comp.h"
#include "spdlog/fmt/fmt.h"
#include "spdlog/fmt/ostr.h"

namespace {

using AlnEngine = std::unique_ptr<spoa::AlignmentEngine>;
inline auto BuildAlnEngine(const usize ref_length, absl::Span<const std::string> alts) -> AlnEngine {
  static constexpr i8 MATCH = 1;
  static constexpr i8 MISMATCH = -9;
  static constexpr i8 OPEN1 = -16;
  static constexpr i8 EXTEND1 = -2;
  static constexpr i8 OPEN2 = -41;
  static constexpr i8 EXTEND2 = -1;
  // asm10 from minimap2 -> https://lh3.github.io/minimap2/minimap2.html -> sequences upto 10% divergence
  // https://curiouscoding.nl/posts/pairwise-alignment -> Convex affine -> min(g1 + (i - 1) * e1, g2 + (i - 1) * e2)
  auto aln = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, MATCH, MISMATCH, OPEN1, EXTEND1, OPEN2, EXTEND2);

  auto max_len = ref_length;
  std::ranges::for_each(alts, [&max_len](const std::string& item) { max_len = std::max(max_len, item.size()); });
  aln->Prealloc(max_len, 4);

  return aln;
}

}  // namespace

namespace lancet::caller {

MsaBuilder::MsaBuilder(RefHaplotype ref_seq, AltHaplotypes alt_seqs, const FsPath& out_gfa_path)
    : mRefRegion(std::move(ref_seq)), mAlternateSeqs(alt_seqs) {
  mResultMsa.reserve(mAlternateSeqs.size() + 1);
  absl::FixedArray<bool> is_reversed(mAlternateSeqs.size() + 1, false);

  spoa::Graph graph;
  const auto aln_engine = BuildAlnEngine(mRefRegion->Length(), mAlternateSeqs);

  AddSequenceToGraph(mRefRegion->SeqData(), mRefRegion->Length(), aln_engine.get(), graph, is_reversed[0]);

  usize curr_idx = 1;
  std::ranges::for_each(mAlternateSeqs, [&is_reversed, &curr_idx, &aln_engine, &graph](const std::string& seq) {
    MsaBuilder::AddSequenceToGraph(seq.c_str(), seq.length(), aln_engine.get(), graph, is_reversed[curr_idx]);
    curr_idx++;
  });

  mResultMsa = graph.GenerateMultipleSequenceAlignment(false);
  if (!out_gfa_path.empty()) {
    WriteFasta(out_gfa_path, mResultMsa);
    WriteGfa(out_gfa_path, graph, is_reversed);
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

void MsaBuilder::WriteGfa(const FsPath& out_path, const spoa::Graph& graph, absl::Span<const bool> is_reversed) {
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
    while (true) {
      path.emplace_back(current_node->id + 1);
      current_node = current_node->Successor(seq_idx);
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (current_node == nullptr) break;
    }

    const auto is_seq_reversed = is_reversed.at(seq_idx);
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (is_seq_reversed) std::reverse(path.begin(), path.end());

    for (usize path_idx = 0; path_idx < path.size(); ++path_idx) {
      fmt::print(out_handle, "{}{}{}", path_idx == 0 ? "" : ",", path[path_idx], is_seq_reversed ? "-" : "+");
    }

    fmt::print(out_handle, "\t*\n");
  }

  out_handle.close();
}

void MsaBuilder::AddSequenceToGraph(const char* sequence, usize length, spoa::AlignmentEngine* aln_engine,
                                    spoa::Graph& graph, bool& added_reverse_complement) {
  i32 fwd_score = 0;
  i32 rev_score = 0;
  const auto rc_seq = RevComp(std::string_view(sequence, length));
  const auto alignment_fwd = aln_engine->Align(sequence, length, graph, &fwd_score);
  const auto alignment_rev = aln_engine->Align(rc_seq, graph, &rev_score);
  const auto is_rev_better = rev_score > fwd_score;
  is_rev_better ? graph.AddAlignment(alignment_rev, rc_seq) : graph.AddAlignment(alignment_fwd, sequence, length);
  added_reverse_complement = is_rev_better;
}

}  // namespace lancet::caller
