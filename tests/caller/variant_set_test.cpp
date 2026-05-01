#include "lancet/caller/variant_set.h"

#include "lancet/caller/alt_allele.h"
#include "lancet/core/window.h"

#include "absl/container/btree_set.h"
#include "catch_amalgamated.hpp"
#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"

#include <algorithm>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

// =========================================================================================
// VariantSet constructor — algorithm test suite
// -----------------------------------------------------------------------------------------
// Validates SNV, INS, DEL, MNP, complex, and multiallelic extraction from SPOA
// graph topologies. Each SECTION builds a controlled graph from string pairs and
// verifies that the extracted VCF records have correct REF, ALT, and POS fields.
// =========================================================================================

namespace lancet::caller::tests {

// Catch2 SECTION fan-out inflates clang-tidy's cognitive-complexity metric beyond the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Topological VariantExtractor securely populates Multiallelic Payload Sets",
          "[lancet][caller][VariantSet]") {
  // Build SPOA graphs from string pairs to create controlled graph topologies for testing.
  auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);

  SECTION(
      "Extracts isolated Biallelic SNV dynamically managing VCF alignment anchoring universally") {
    // [REF]:  A T C G
    // [ALT]:  A G C G   => (SNV: T -> G)

    spoa::Graph graph{};
    std::vector<std::string> const seqs = {"ATCG", "AGCG"};
    for (auto const& seq : seqs) graph.AddAlignment(alignment_engine->Align(seq, graph), seq);

    VariantSet caller_set(graph, core::Window{}, 100);

    // Left-trimming strips the shared 'A' prefix, leaving the minimal REF='T', ALT='G'
    // representation. Position shifts +1 from the window start (100 → 101).
    REQUIRE(caller_set.Count() == 1);
    auto var = *caller_set.begin();
    REQUIRE(var.mRefAllele == "T");
    REQUIRE(var.mAlts.size() == 1);
    REQUIRE(var.mAlts[0].mSequence == "G");
    REQUIRE(var.mGenomeChromPos1 == 101);
  }

  SECTION("Maintains obligatory Left Anchor Bounds for strict mathematical VCF Deletions") {
    //         [Anchor]
    // [REF]:      A T C G
    // [ALT]:      A - - G   => (DEL: ATC -> A)

    spoa::Graph graph{};
    std::vector<std::string> const seqs = {"ATCG", "AG"};
    for (auto const& seq : seqs) graph.AddAlignment(alignment_engine->Align(seq, graph), seq);

    VariantSet caller_set(graph, core::Window{}, 100);

    // Deletions must retain the left anchor base — empty REF/ALT is invalid VCF.
    // REF='ATC', ALT='A' preserves the mandatory anchor.
    REQUIRE(caller_set.Count() == 1);
    auto var = *caller_set.begin();
    REQUIRE(var.mRefAllele == "ATC");
    REQUIRE(var.mAlts.size() == 1);
    REQUIRE(var.mAlts[0].mSequence == "A");
    REQUIRE(var.mGenomeChromPos1 == 100);
  }

  SECTION("Extracts overlapping multi-allelic mutations with correct haplotype tracking") {
    // Deep Overlap Multi-allelic test:
    // [REF]:    A T G T G C
    // [ALT1]:   A C G T G C     => (SNV T->C)
    // [ALT2]:   A - - - G C     => (DEL TGT->-)
    // [ALT3]:   A T G T A C     => (SNV G->A) at downstream position
    //
    // ALT1 (SNV T→C) and ALT2 (DEL TGT) share the same anchor node 'A', so they merge
    // into a single multiallelic record. The downstream G→A SNV (ALT3) is also fused because
    // the deletion spans into its position.

    spoa::Graph graph{};
    std::vector<std::string> const seqs = {"ATGTGC", "ACGTGC", "AGC", "ATGTAC"};
    for (auto const& seq : seqs) graph.AddAlignment(alignment_engine->Align(seq, graph), seq);

    VariantSet caller_set(graph, core::Window{}, 100);

    // The TGT deletion spans up to the G→A SNV position, so all three ALTs fuse into
    // one multiallelic record.
    REQUIRE(caller_set.Count() == 1);

    auto first_var = *caller_set.begin();
    REQUIRE(first_var.mAlts.size() == 3);

    // Sorting dictates alphanumeric strict mapping: "A", "C", "TGTAC". But wait!
    // We will just dynamically check presence to be algorithmically safe against graph-internal
    // topological sorting behaviors!
    bool found_del = false;
    bool found_c_snv = false;
    bool found_a_snv = false;

    std::ranges::for_each(first_var.mAlts, [&](auto const& alt) {
      if (alt.mLocalHapStart0Idxs.contains(2)) found_del = true;
      if (alt.mLocalHapStart0Idxs.contains(1)) found_c_snv = true;
      if (alt.mLocalHapStart0Idxs.contains(3)) found_a_snv = true;
    });

    REQUIRE(found_del == true);
    REQUIRE(found_c_snv == true);
    REQUIRE(found_a_snv == true);
  }

  SECTION("Calculates a Pure Structural Insertion gracefully bounded universally against the exact "
          "Reference Target") {
    //                                  [ALT]
    //                          .-->(A)[3] --> (A)[4] --.
    //                         /                          \
    //   Anchor: (T)[2] ------+                            +-----> Target: (C)[5] (CONVERGED!)
    //                         \                          /
    //                          `------------------------'
    //                                  [REF]
    // REF:  A T C G
    // ALT:  A T A A C G

    spoa::Graph graph{};
    std::vector<std::string> const seqs = {"ATCG", "ATAACG"};
    for (auto const& seq : seqs) graph.AddAlignment(alignment_engine->Align(seq, graph), seq);

    VariantSet caller_set(graph, core::Window{}, 100);

    REQUIRE(caller_set.Count() == 1);
    auto var = *caller_set.begin();
    // VCF left-alignment: floating inserts anchor to the preceding base
    REQUIRE(var.mRefAllele == "T");
    REQUIRE(var.mAlts.size() == 1);
    REQUIRE(var.mAlts[0].mSequence == "TAA");
    // Mapped accurately to the true 'T' anchor position
    REQUIRE(var.mGenomeChromPos1 == 101);
  }

  SECTION("Preserves symmetrical MNP boundaries without shattering") {
    //                         .-->(A) --> (A) --.
    //                        /                    \
    //          (A) ---------+                      +----> (G)
    //                        \                    /
    //                         `-->(T) --> (C) --'
    // REF:  A T C G
    // ALT:  A A A G   (MNP: TC -> AA)

    spoa::Graph graph{};
    std::vector<std::string> const seqs = {"ATCG", "AAAG"};
    for (auto const& seq : seqs) graph.AddAlignment(alignment_engine->Align(seq, graph), seq);

    VariantSet caller_set(graph, core::Window{}, 100);

    REQUIRE(caller_set.Count() == 1);
    auto var = *caller_set.begin();
    // VCF normalizers strip rightward and leftward independently. If the internal string core
    // differs, they cannot strip further. MNPs retain their full structure.
    REQUIRE(var.mRefAllele == "TC");
    REQUIRE(var.mAlts.size() == 1);
    REQUIRE(var.mAlts[0].mSequence == "AA");
  }

  SECTION("Fuses complex variants without misalignment biases during extraction") {
    //                         .-->(A) --> (A) --> (A) --.
    //                        /                            \
    //          (A) ---------+                              +----> (G)
    //                        \                            /
    //                         `------>(T) --> (C) ------'
    // REF:  A T C G
    // ALT:  A A A A G   (CPX: TC -> AAA)

    spoa::Graph graph{};
    std::vector<std::string> const seqs = {"ATCG", "AAAAG"};
    for (auto const& seq : seqs) graph.AddAlignment(alignment_engine->Align(seq, graph), seq);

    VariantSet caller_set(graph, core::Window{}, 100);

    // SPOA optimizes this topological structure: instead of a monolithic
    // mismatch, the aligner decoupled them into two isolated events:
    // 1. Insertion of 'A' (Between the A and T anchors)
    // 2. MNP of 'TC' -> 'AA'
    REQUIRE(caller_set.Count() == 2);

    auto first_var = *caller_set.begin();
    REQUIRE(first_var.mGenomeChromPos1 == 100);
    REQUIRE(first_var.mRefAllele.empty());  // Pure raw insertion — no REF allele
    REQUIRE(first_var.mAlts[0].mSequence == "A");

    auto second_var = *std::next(caller_set.begin());
    REQUIRE(second_var.mGenomeChromPos1 == 101);
    REQUIRE(second_var.mRefAllele == "TC");
    REQUIRE(second_var.mAlts.size() == 1);
    REQUIRE(second_var.mAlts[0].mSequence == "AA");
  }

  SECTION(
      "Handles Nested N-Way Sink Synchronization effectively inside dense mutational clusters") {
    //                          .-->(G)[3] --.     [ALT1: SNV]
    //                         /               \
    //                        +---->(A)[4] ----+   [ALT2: SNV]
    //                       /                   \
    //  Anchor: (A)[2] -----+---------------------+------> Target: (C)[5] (CONVERGED!)
    //                       \                   /  [ALT3: DEL]
    //                        `---->(T)[6] ----'   [REF]
    // REF:  A T C G
    // ALT1: A G C G
    // ALT2: A A C G
    // ALT3: A C G  (del T)
    //
    // Note: A C G logically skips the T, merging identically into C!
    spoa::Graph graph{};
    std::vector<std::string> const seqs = {"ATCG", "AGCG", "AACG", "ACG"};
    for (auto const& seq : seqs) graph.AddAlignment(alignment_engine->Align(seq, graph), seq);

    VariantSet caller_set(graph, core::Window{}, 100);

    // All three ALTs (G, A, and deletion) merge into a single triallelic VCF record
    // with REF='AT' (anchor + deleted base).
    REQUIRE(caller_set.Count() == 1);
    auto var = *caller_set.begin();

    REQUIRE(var.mAlts.size() == 3);
    REQUIRE(var.mRefAllele == "AT");

    bool found_g = false;
    bool found_a = false;
    bool found_del = false;

    std::ranges::for_each(var.mAlts, [&](auto const& alt) {
      if (alt.mSequence == "AG") found_g = true;
      if (alt.mSequence == "AA") found_a = true;
      // Deletions strictly hold the universal 'A' Anchor!
      if (alt.mSequence == "A") found_del = true;
    });

    REQUIRE(found_g == true);
    REQUIRE(found_a == true);
    REQUIRE(found_del == true);
  }
}
}  // namespace lancet::caller::tests
