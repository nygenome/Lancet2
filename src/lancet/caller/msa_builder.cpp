#include "lancet/caller/msa_builder.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <ios>
#include <iterator>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "spdlog/fmt/bundled/ostream.h"
#include "spoa/alignment_engine.hpp"

namespace {

/*
 * ============================================================================
 * SPOA MSA Parameter Rationale for Lancet2 Variant Extraction
 * ============================================================================
 * Values: Match: 2, Mismatch: -4, Gap1: -4,-2, Gap2: -24,-1
 * 
 * Unlike minimap2's `asm5` preset (which aggressively splits contigs at major 
 * divergences for whole-genome synteny filtering), these parameters are tuned 
 * to force end-to-end global alignment within a specific micro-assembly window 
 * to capture dense somatic mutations and large Insertions/Deletions.
 * 
 * 1. Why Convex (Dual-Affine) vs. Affine or Linear Scoring:
 *    - Linear Scoring applies a flat penalty per gap base, which is biologically 
 *      inaccurate (one 50bp deletion is one biological event, not fifty 1bp 
 *      independent events).
 *    - Single Affine Scoring forces a compromise: tune for small variants (strict 
 *      extension) and you penalize/clip large insertions/deletions; tune for large 
 *      insertions/deletions (loose extension) and sequencer noise creates messy, 
 *      spurious small gaps.
 *    - Convex (Dual-Affine) Scoring solves this by taking the minimum of two 
 *      intersecting models. It is strict for short gaps to suppress sequencer 
 *      noise, but switches to an incredibly cheap extension penalty for large 
 *      biological gaps.
 * 
 * 2. Mismatch Tolerance (Multi-Nucleotide Variants / MNVs): 
 *    asm5 uses a +1 match / -19 mismatch, which shatters alignments at dense 
 *    mutation clusters. We use +2 / -4, meaning only 2 matching bases are needed 
 *    to offset a SNP. This keeps the MSA globally intact through complex variants.
 * 
 * 3. Micro-Indel Sensitivity (Convex Model 1: -4, -2): 
 *    asm5's -39 gap open penalty prevents small indels, forcing them to misalign 
 *    as false-positive SNPs. Our -4 open / -2 extend penalty allows true small 
 *    biological indels to open naturally while still applying enough friction 
 *    to prevent 1bp sequencing errors (e.g., homopolymer stutters) from opening gaps.
 * 
 * 4. Large Insertion/Deletion Continuity (Convex Model 2: -24, -1): 
 *    asm5's -81 penalty for large gaps will soft-clip contigs right at an insertion/deletion 
 *    breakpoint. Our parameters mathematically intersect at 20bp (4 + 2L = 24 + 1L). 
 *    For gaps > 20bp, the algorithm switches to Model 2 where the extension cost 
 *    drops to -1. This "cheap extension" forces the DP matrix into mapping massive 
 *    insertions/deletions as single, contiguous blocks in the MSA rather than 
 *    dropping the alignment.
 * 
 *    https://curiouscoding.nl/posts/pairwise-alignment -> 
 *    – Convex dual affine gap scoring -> min(g1+(i-1)*e1, g2+(i-1)*e2)
 */
using AlnEngine = std::unique_ptr<spoa::AlignmentEngine>;
inline auto BuildAlnEngine(absl::Span<const std::string> seqs) -> AlnEngine {
  static constexpr i8 MATCH = 2;
  static constexpr i8 MISMATCH = -4;
  static constexpr i8 OPEN1 = -4;
  static constexpr i8 EXTEND1 = -2;
  static constexpr i8 OPEN2 = -24;
  static constexpr i8 EXTEND2 = -1;
  auto aln = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, MATCH, MISMATCH, OPEN1, EXTEND1, OPEN2, EXTEND2);
  const auto* longest_seq = std::ranges::max_element(seqs, std::less<>(), &std::string::size);
  aln->Prealloc(longest_seq->length(), 16);
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

    std::vector<u32> path;
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
