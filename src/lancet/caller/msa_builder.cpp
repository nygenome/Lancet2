#include "lancet/caller/msa_builder.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <ios>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>

#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "spdlog/fmt/bundled/ostream.h"
#include "spoa/alignment_engine.hpp"


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

void MsaBuilder::UpdateSpoaState(absl::Span<const std::string> sequences) {
  mGraph.Clear();
  std::ranges::for_each(sequences, [this](const std::string& haplotype) {
    const auto alignment = mEngine->Align(haplotype, mGraph);
    mGraph.AddAlignment(alignment, haplotype);
  });
}

void MsaBuilder::SerializeGraph(const FsPath& out_gfa_path) {
  if (out_gfa_path.empty()) return;
  const auto msa_alns = mGraph.GenerateMultipleSequenceAlignment(false);
  WriteFasta(out_gfa_path, msa_alns);
  WriteGfa(out_gfa_path);
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

void MsaBuilder::WriteGfa(const FsPath& out_path) const {
  // https://github.com/rvaser/spoa/pull/36
  // See PR for how to normalize & process the output GFA
  std::ofstream out_handle(out_path, std::ios::trunc);
  fmt::print(out_handle, "H\tVN:Z:1.0\n");

  // --- 1. Write Segments (Nodes) and Links (Edges) ---
  for (const std::unique_ptr<spoa::Graph::Node>& node : mGraph.nodes()) {
    const auto node_id = node->id + 1;
    const auto node_seq = static_cast<char>(mGraph.decoder(node->code));
    fmt::print(out_handle, "S\t{}\t{}\n", node_id, node_seq);
    
    for (const spoa::Graph::Edge* edge : node->outedges) {
      const auto dest_id = edge->head->id + 1;
      fmt::print(out_handle, "L\t{}\t+\t{}\t+\t0M\n", node_id, dest_id);
    }
  }

  // --- 2. Write Walks (Sequence Paths) ---
  for (u32 seq_idx = 0; seq_idx < mGraph.sequences().size(); ++seq_idx) {
    const auto hap_prefix = seq_idx == 0 ? "ref" : "hap";
    fmt::print(out_handle, "P\t{}{}\t", hap_prefix, seq_idx);

    // Iteratively trace node Successors to construct the sequence path organically
    const spoa::Graph::Node* curr_node = mGraph.sequences()[seq_idx];
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
