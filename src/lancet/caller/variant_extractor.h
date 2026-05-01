#ifndef SRC_LANCET_CALLER_VARIANT_EXTRACTOR_H_
#define SRC_LANCET_CALLER_VARIANT_EXTRACTOR_H_

#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_bubble.h"
#include "lancet/core/window.h"

#include "absl/container/btree_set.h"
#include "absl/types/span.h"
#include "spoa/graph.hpp"

#include <limits>
#include <string>
#include <vector>

namespace lancet::caller {

constexpr usize REF_HAP_IDX = 0;

// =========================================================================================
// TOPOLOGICAL BUBBLE EXTRACTION — How Variants Are Discovered in the POA Graph
// =========================================================================================
// Rather than relying on 2D string matrices (MSA column diffing)—which scale poorly—or
// pairwise mapping—which fractures overlapping multiallelic loci—this algorithm tracks
// pointers traversing the SPOA directed acyclic graph (DAG).
//
// -> 1. BIALLELIC SNVs (Paths diverge at a single node, then reconverge):
//
//                                [REF]
//                          .--> (T)[3] --.
//                         /               \
//   Anchor: (A)[2] ------+                 +-----> Target: (G)[5] (CONVERGED!)
//                         \               /
//                          `--> (C)[4] --'
//                                [ALT]
//
// -> 2. DELETIONS (ALT path skips a stretch of REF nodes via a direct edge):
//
//                                            [REF]
//                          .--> (T)[3] --> (C)[4] --> (G)[5] --.
//                         /                                     \
//   Anchor: (A)[2] ------+                                       +-----> Target: (T)[6]
//   (CONVERGED!)
//                         \                                     /
//                          `-----------------------------------'
//                                            [ALT]
//
// -> 3. MULTIALLELIC COMPLEXES (Multiple ALT paths diverge simultaneously):
//       Resolving these requires advancing all path pointers in topological order.
//
//                               [ALT 1]
//                          .--> (C)[3] --> (A)[4] ---------.
//                         /                                 \
//   Anchor: (T)[2] ------+----> (T)[5] --> (G)[6] --> (C)[7] +-----> Target: (A)[10] (CONVERGED!)
//                         \     [REF]                       /
//                          `--> (G)[8] --> (A)[9] ---------'
//                               [ALT 2]
//
// HOW THE SWEEP EXTRACTOR ALGORITHM WORKS ("Greedy Sink"):
// 1. Maintain an array of active pointers — one per haplotype — tracking the
//    current node each path occupies in the DAG.
// 2. When pointers disagree (point to different nodes), a "Bubble" has opened.
// 3. Each pointer has a `Rank` — its topologically sorted 5'-to-3' position
//    in the graph. Rank is the linearized left-to-right ordering.
// 4. Repeatedly advance only the pointer(s) with the LOWEST rank, appending
//    their decoded bases to per-path string buffers.
// 5. This causes lagging paths to catch up until ALL pointers converge on
//    the same Rank — the bubble's convergence node.
// 6. Collect the buffered sequences, pass them through VCF multi-allelic
//    trimming (NormalizeVcfParsimony), and emit the RawVariant.
// ================================================================================================
class VariantExtractor {
 public:
  VariantExtractor(spoa::Graph const& graph, core::Window const& win, usize anchor_start);

  // ===================================================================================================
  // VARIANT EXTRACTOR FLOWCHART
  // ===================================================================================================
  // Pushes active pointers 5'-to-3' through the POA graph. When pathways diverge
  // (a variant occurs), the extractor halts uniform sweeping and triggers a Bubble.
  //
  //                              (C)  (Prior Match Node — all pointers converged here)
  //                            /     \
  //                  (ALT)    /       \   (REF)
  //                 (C)->(T) .         . (T)->(G)
  //                           \       /
  //                            \     /
  //                              (G)      (Target Convergence Node — paths reunite here)
  //
  // Private helper call sequence:
  // 1. `InitializeBubbleAnchor`    : Prepends the last matched base (C) as the VCF anchor.
  // 2. `SinkPointers`              : The hot loop. Finds the lowest-rank pointer,
  //    |-- `FindLowestActiveRank`  :   advances it, and repeats until all pointers
  //    |-- `ConsumePathsAtRank`    :   converge on (G).
  // 3. `CreateNormalizedBubble`    : Groups per-path strings, runs VCF parsimony trimming.
  // 4. `AssembleMultiallelicVariant`: Classifies each ALT and emits a RawVariant.
  // ===================================================================================================
  void SearchAndExtractTo(absl::btree_set<RawVariant>& out_variants);

 private:
  // VariantExtractor is a transient extraction context bound to a single (graph, window) pair —
  // const-ref members are deliberate and the class is non-copyable by intent.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::vector<u32> mNodeToRank;                       // 8B (24B)
  std::vector<spoa::Graph::Node const*> mActivePtrs;  // 8B (24B)
  std::vector<usize> mCurrentHapPos;                  // 8B (24B)
  spoa::Graph const& mGraph;                          // 8B
  core::Window const& mWin;                           // 8B
  usize mRefAnchorStart;                              // 8B
  usize mCurrentRefPos;                               // 8B
  usize mNumSeqs;                                     // 8B

  // Last converged VCF anchor node (prepended to bubble sequences)
  spoa::Graph::Node const* mPrevMatchNode = nullptr;  // 8B
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

  // True when all haplotype pointers point to the same DAG node.
  [[nodiscard]] auto AreAllPathsConverged() const -> bool;

  // Step all converged pointers to their next successor node.
  void AdvanceConvergedPaths();

  // Detect a single bubble divergence and emit its variant(s).
  void EatTopologicalBubble(absl::btree_set<RawVariant>& out_variants);

  // Prepend the last matched base as VCF anchor; returns bubble start position.
  auto InitializeBubbleAnchor(absl::Span<std::string> raw_alleles,
                              std::vector<usize>& out_hap_starts) -> usize;

  // Advance lagging pointers until all paths reconverge.
  void SinkPointers(absl::Span<std::string> raw_alleles);

  // Return the minimum topological rank among active (non-null) pointers.
  [[nodiscard]] auto FindLowestActiveRank() const -> u32;

  // Advance pointers at the given rank, appending decoded bases to allele strings.
  void ConsumePathsAtRank(u32 target_rank, absl::Span<std::string> raw_alleles);

  // Group per-path sequences and apply VCF parsimony trimming.
  [[nodiscard]] auto CreateNormalizedBubble(usize genome_start_pos,
                                            std::vector<std::string> raw_alleles) const
      -> VariantBubble;

  // Classify each ALT allele and build the final multiallelic RawVariant.
  auto AssembleMultiallelicVariant(VariantBubble bubble) -> RawVariant;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_EXTRACTOR_H_
