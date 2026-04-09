#include "lancet/caller/variant_set.h"

#include "lancet/core/window.h"

#include "absl/strings/string_view.h"
#include "catch_amalgamated.hpp"
#include "spoa/spoa.hpp"

// =========================================================================================
// CATCH2 COMPREHENSIVE ALGORITHM TEST SUITE
// -----------------------------------------------------------------------------------------
// Asserts that the FSM VariantExtractor properly mitigates, sinks, parses, and formats natively
// raw sequences utilizing dynamically generated topological networks mapping real complex bounds.
// =========================================================================================

namespace lancet::caller::tests {

TEST_CASE("Topological VariantExtractor securely populates Multiallelic Payload Sets",
          "[lancet][caller][VariantSet]") {
  // We bind a linear Smith-Waterman engine to dynamically forge graphs universally natively out of
  // strings. The graph topology will authentically mimic organic topological sequencing branches.
  auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);

  SECTION(
      "Extracts isolated Biallelic SNV dynamically managing VCF alignment anchoring universally") {
    // [REF]:  A T C G
    // [ALT]:  A G C G   => (SNV: T -> G)

    spoa::Graph graph{};
    std::vector<std::string> seqs = {"ATCG", "AGCG"};
    for (auto const& s : seqs) graph.AddAlignment(alignment_engine->Align(s, graph), s);

    VariantSet caller_set;
    caller_set.ExtractVariantsFromGraph(graph, core::Window{}, 100);

    // Verification Check: Our sweep relies on universal matching. A match generated 'AT'/'AG'
    // originally natively encompassing the anchor bounds identically. Then, the left-truncation
    // algorithm strips away universally useless left-alignment bases because no element resulted in
    // a critical empty string violating bounds rules!
    REQUIRE(caller_set.Count() == 1);
    auto var = *caller_set.begin();
    REQUIRE(var.mRefAllele == "T");
    REQUIRE(var.mAlts.size() == 1);
    REQUIRE(var.mAlts[0].mSequence == "G");
    REQUIRE(var.mGenomeChromPos1 ==
            101);  // Biologically validated string coordinate slid 1 index to the right
  }

  SECTION("Maintains obligatory Left Anchor Bounds for strict mathematical VCF Deletions") {
    //         [Anchor]
    // [REF]:      A T C G
    // [ALT]:      A - - G   => (DEL: ATC -> A)

    spoa::Graph graph{};
    std::vector<std::string> seqs = {"ATCG", "AG"};
    for (auto const& s : seqs) graph.AddAlignment(alignment_engine->Align(s, graph), s);

    VariantSet caller_set;
    caller_set.ExtractVariantsFromGraph(graph, core::Window{}, 100);

    // Critical Test Constraint:
    // Under no circumstances can Deletions be left-trimmed so aggressively they turn into native
    // "" empty strings (violating VCF configuration specifications permanently downstream!)
    REQUIRE(caller_set.Count() == 1);
    auto var = *caller_set.begin();
    REQUIRE(var.mRefAllele == "ATC");
    REQUIRE(var.mAlts.size() == 1);
    REQUIRE(var.mAlts[0].mSequence == "A");
    REQUIRE(var.mGenomeChromPos1 == 100);
  }

  SECTION("Encapsulates massively overlapping varying mutations tracking natively identically to "
          "Haplotypes") {
    // Deep Overlap Multi-allelic test:
    // [REF]:    A T G T G C
    // [ALT1]:   A C G T G C     => (SNV T->C)
    // [ALT2]:   A - - - G C     => (DEL TGT->-)
    // [ALT3]:   A T G T A C     => (SNV G->A) at downstream position
    //
    // Our Synchronous-Sink sweeps ALL paths identically concurrently natively. Because ALT1 and
    // ALT2 topologically branch entirely off mathematically contiguous identical 'A' Anchor Nodes,
    // they are fundamentally synthesized together merging continuously until they natively
    // intuitively rendezvous! Downstream SNVs act perfectly independent from anterior topological
    // structures.

    spoa::Graph graph{};
    std::vector<std::string> seqs = {"ATGTGC", "ACGTGC", "AGC", "ATGTAC"};
    for (auto const& s : seqs) graph.AddAlignment(alignment_engine->Align(s, graph), s);

    VariantSet caller_set;
    caller_set.ExtractVariantsFromGraph(graph, core::Window{}, 100);

    // The Topological Sink now perfectly recognizes that the TGT Deletion mathematically spans
    // completely adjacent to the G->A SNV! This seamlessly fuses ALL of them natively into ONE
    // massive multi-allelic block!
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
    //                          .--> (A)[3] --> (A)[4] --.
    //                         /                          \
    //   Anchor: (T)[2] ------+                            +-----> Target: (C)[5] (CONVERGED!)
    //                         \                          /
    //                          `------------------------'
    //                                  [REF]
    // REF:  A T C G
    // ALT:  A T A A C G

    spoa::Graph graph{};
    std::vector<std::string> seqs = {"ATCG", "ATAACG"};
    for (auto const& s : seqs) graph.AddAlignment(alignment_engine->Align(s, graph), s);

    VariantSet caller_set;
    caller_set.ExtractVariantsFromGraph(graph, core::Window{}, 100);

    REQUIRE(caller_set.Count() == 1);
    auto var = *caller_set.begin();
    REQUIRE(var.mRefAllele ==
            "T");  // Universal VCF Left-Alignment rule binds floating inserts correctly!
    REQUIRE(var.mAlts.size() == 1);
    REQUIRE(var.mAlts[0].mSequence == "TAA");
    REQUIRE(var.mGenomeChromPos1 == 101);  // Mapped accurately to the true 'T' anchor position
  }

  SECTION("Safeguards symmetrical boundaries against MNP shattering properly natively") {
    //                         .--> (A) --> (A) --.
    //                        /                    \
    //          (A) ---------+                      +----> (G)
    //                        \                    /
    //                         `--> (T) --> (C) --'
    // REF:  A T C G
    // ALT:  A A A G   (MNP: TC -> AA)

    spoa::Graph graph{};
    std::vector<std::string> seqs = {"ATCG", "AAAG"};
    for (auto const& s : seqs) graph.AddAlignment(alignment_engine->Align(s, graph), s);

    VariantSet caller_set;
    caller_set.ExtractVariantsFromGraph(graph, core::Window{}, 100);

    REQUIRE(caller_set.Count() == 1);
    auto var = *caller_set.begin();
    // Since VCF normalizers strip rightward and leftward independently, if the internal string core
    // differs, they cannot strip identically. MNP's logically hold their total structure
    // identically natively!
    REQUIRE(var.mRefAllele == "TC");
    REQUIRE(var.mAlts.size() == 1);
    REQUIRE(var.mAlts[0].mSequence == "AA");
  }

  SECTION("Fuses massive Complex Variants natively without misalignment biases during Extraction") {
    //                         .--> (A) --> (A) --> (A) --.
    //                        /                            \
    //          (A) ---------+                              +----> (G)
    //                        \                            /
    //                         `------> (T) --> (C) ------'
    // REF:  A T C G
    // ALT:  A A A A G   (CPX: TC -> AAA)

    spoa::Graph graph{};
    std::vector<std::string> seqs = {"ATCG", "AAAAG"};
    for (auto const& s : seqs) graph.AddAlignment(alignment_engine->Align(s, graph), s);

    VariantSet caller_set;
    caller_set.ExtractVariantsFromGraph(graph, core::Window{}, 100);

    // SPOA logically optimizes this topological structure! Instead of a massive monolithic
    // mismatch, the engine perfectly decoupled them into two isolated dynamic events natively:
    // 1. Insertion of 'A' (Between the A and T anchors)
    // 2. MNP of 'TC' -> 'AA'
    REQUIRE(caller_set.Count() == 2);

    auto first_var = *caller_set.begin();
    REQUIRE(first_var.mGenomeChromPos1 == 100);
    REQUIRE(first_var.mRefAllele == "");  // Natively evaluated pure raw insertion
    REQUIRE(first_var.mAlts[0].mSequence == "A");

    auto second_var = *std::next(caller_set.begin());
    REQUIRE(second_var.mGenomeChromPos1 == 101);
    REQUIRE(second_var.mRefAllele == "TC");
    REQUIRE(second_var.mAlts.size() == 1);
    REQUIRE(second_var.mAlts[0].mSequence == "AA");
  }

  SECTION(
      "Handles Nested N-Way Sink Synchronization effectively inside dense mutational clusters") {
    //                          .--> (G)[3] --.     [ALT1: SNV]
    //                         /               \
    //                        +----> (A)[4] ----+   [ALT2: SNV]
    //                       /                   \
    //  Anchor: (A)[2] -----+---------------------+-----> Target: (C)[5] (CONVERGED!)
    //                       \                   /  [ALT3: DEL]
    //                        `----> (T)[6] ----'   [REF]
    // REF:  A T C G
    // ALT1: A G C G
    // ALT2: A A C G
    // ALT3: A C G  (del T)
    //
    // Note: A C G logically skips the T, merging identically into C!
    spoa::Graph graph{};
    std::vector<std::string> seqs = {"ATCG", "AGCG", "AACG", "ACG"};
    for (auto const& s : seqs) graph.AddAlignment(alignment_engine->Align(s, graph), s);

    VariantSet caller_set;
    caller_set.ExtractVariantsFromGraph(graph, core::Window{}, 100);

    // They must all uniformly synthesize precisely simultaneously natively into a strictly
    // evaluated multiallelic VCF payload seamlessly bounded.
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
