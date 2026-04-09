#include "lancet/caller/msa_builder.h"

#include "lancet/base/types.h"

#include "absl/types/span.h"
#include "spdlog/fmt/bundled/ostream.h"
#include "spoa/alignment_engine.hpp"

#include <algorithm>
#include <fstream>
#include <ios>
#include <memory>
#include <string>
#include <vector>

/*
 * NOTE: Both the AlignmentEngine and the spoa::Graph are instantiated once per
 * worker thread within the caller::MsaBuilder object in VariantBuilder. They are
 * passed here by reference. The engine maintains the scoring params and preallocated
 * DP matrix buffer. The graph, which stores the topological biological state,
 * is explicitly Clear()ed per window component. This eliminates the massive OS heap
 * fragmentation penalty from continuous graph initialization/destruction, ensuring
 * strict zero-allocation maximum execution speed across millions of genomic windows.
 */

namespace lancet::caller {

void MsaBuilder::UpdateSpoaState(absl::Span<std::string const> sequences) {
  mGraph.Clear();
  std::ranges::for_each(sequences, [this](std::string const& haplotype) -> void {
    auto const alignment = mEngine->Align(haplotype, mGraph);
    mGraph.AddAlignment(alignment, haplotype);
  });
}

void MsaBuilder::SerializeGraph(FsPath const& out_gfa_path) {
  if (out_gfa_path.empty()) {
    return;
  }
  auto const msa_alns = mGraph.GenerateMultipleSequenceAlignment(false);
  WriteFasta(out_gfa_path, msa_alns);
  WriteGfa(out_gfa_path);
}

void MsaBuilder::WriteFasta(FsPath const& gfa_path, absl::Span<std::string const> msa_alns) {
  auto fa_path = gfa_path.parent_path() / gfa_path.stem();
  fa_path += ".fasta";
  std::ofstream out_handle(fa_path);
  for (usize idx = 0; idx < msa_alns.size(); ++idx) {
    fmt::print(out_handle, ">{}{}\n{}\n", idx == 0 ? "ref" : "hap", idx, msa_alns[idx]);
  }
  out_handle.close();
}

void MsaBuilder::WriteGfa(FsPath const& out_path) const {
  // https://github.com/rvaser/spoa/pull/36
  // See PR for how to normalize & process the output GFA
  std::ofstream out_handle(out_path, std::ios::trunc);
  fmt::print(out_handle, "H\tVN:Z:1.0\n");

  // --- 1. Write Segments (Nodes) and Links (Edges) ---
  for (std::unique_ptr<spoa::Graph::Node> const& node : mGraph.nodes()) {
    auto const node_id = node->id + 1;
    auto const node_seq = static_cast<char>(mGraph.decoder(node->code));
    fmt::print(out_handle, "S\t{}\t{}\n", node_id, node_seq);

    for (spoa::Graph::Edge const* edge : node->outedges) {
      auto const dest_id = edge->head->id + 1;
      fmt::print(out_handle, "L\t{}\t+\t{}\t+\t0M\n", node_id, dest_id);
    }
  }

  // --- 2. Write Walks (Sequence Paths) ---
  for (u32 seq_idx = 0; seq_idx < mGraph.sequences().size(); ++seq_idx) {
    auto const* const hap_prefix = seq_idx == 0 ? "ref" : "hap";
    fmt::print(out_handle, "P\t{}{}\t", hap_prefix, seq_idx);

    // Iteratively trace node Successors to construct the sequence path organically
    spoa::Graph::Node const* curr_node = mGraph.sequences()[seq_idx];
    bool is_first = true;
    while (curr_node != nullptr) {
      fmt::print(out_handle, "{}{}+", is_first ? "" : ",", curr_node->id + 1);
      is_first = false;
      curr_node = curr_node->Successor(seq_idx);
    }

    fmt::print(out_handle, "\t*\n");
  }
}

}  // namespace lancet::caller
