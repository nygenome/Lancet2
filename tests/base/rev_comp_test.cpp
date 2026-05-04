#include "lancet/base/rev_comp.h"

#include "catch_amalgamated.hpp"

#include <initializer_list>
#include <string>
#include <string_view>

namespace lancet::base::tests {

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  RevComp(char) — single-base overload                                    ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("RevComp(char) maps the four canonical bases to their complements",
          "[lancet][base][RevComp]") {
  // The Watson-Crick pairings are A↔T and C↔G. The complement table is built
  // at compile time, so this test pins the ground-truth mapping that every
  // higher-level test (graph traversal, k-mer canonicalization) silently
  // depends on.
  CHECK(RevComp('A') == 'T');
  CHECK(RevComp('T') == 'A');
  CHECK(RevComp('C') == 'G');
  CHECK(RevComp('G') == 'C');
}

TEST_CASE("RevComp(char) preserves case for canonical bases", "[lancet][base][RevComp]") {
  // Some upstream readers emit lowercase soft-masked bases (e.g. RepeatMasker
  // output). Preserving case lets downstream filters distinguish soft-masked
  // alignments from confidently-mapped bases without re-masking.
  CHECK(RevComp('a') == 't');
  CHECK(RevComp('t') == 'a');
  CHECK(RevComp('c') == 'g');
  CHECK(RevComp('g') == 'c');
}

TEST_CASE("RevComp(char) maps N (and n) to itself", "[lancet][base][RevComp]") {
  // N is the IUPAC code for "any base". Its complement is also "any base", so
  // the table maps N→N (and n→n) explicitly rather than falling through to
  // the default-N path used for non-DNA characters.
  CHECK(RevComp('N') == 'N');
  CHECK(RevComp('n') == 'n');
}

TEST_CASE("RevComp(char) maps non-DNA characters to N", "[lancet][base][RevComp]") {
  // Any byte that isn't A/C/G/T/N (case-insensitive) is replaced with N. This
  // turns IUPAC ambiguity codes (R/Y/S/W/K/M/B/D/H/V) and stray whitespace
  // into a defined sentinel rather than passing the original byte through —
  // a regression that returned the original char would silently break
  // strand-symmetry assumptions downstream.
  CHECK(RevComp('R') == 'N');  // IUPAC purine ambiguity
  CHECK(RevComp('Y') == 'N');  // IUPAC pyrimidine ambiguity
  CHECK(RevComp('X') == 'N');  // not a DNA code at all
  CHECK(RevComp(' ') == 'N');
  CHECK(RevComp('\0') == 'N');
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  RevComp(string_view) — sequence overload                                ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("RevComp(string_view) reverse-complements canonical DNA", "[lancet][base][RevComp]") {
  // Hand-checked: ACGT → reverse → TGCA → complement → ACGT. Identity on this
  // input is the standard sanity check for the reverse-then-complement order.
  // The second case (ATCG → CGAT) actually changes; ACGT is a palindrome.
  CHECK(RevComp(std::string_view("ACGT")) == "ACGT");
  CHECK(RevComp(std::string_view("ATCG")) == "CGAT");
  CHECK(RevComp(std::string_view("AAAA")) == "TTTT");
  CHECK(RevComp(std::string_view("GATTACA")) == "TGTAATC");
}

TEST_CASE("RevComp(string_view) returns empty on empty input", "[lancet][base][RevComp]") {
  // Empty input is a legitimate edge case (e.g. a zero-length variant
  // window). The implementation must not deref a null pointer or under-flow
  // the size_t arithmetic.
  CHECK(RevComp(std::string_view("")).empty());
}

TEST_CASE("RevComp(string_view) handles single-base input", "[lancet][base][RevComp]") {
  // Smallest non-empty sequence; the reverse step is a no-op so this case
  // pins the complement-only path.
  CHECK(RevComp(std::string_view("A")) == "T");
  CHECK(RevComp(std::string_view("c")) == "g");
}

TEST_CASE("RevComp(string_view) maps an all-N input to all-N", "[lancet][base][RevComp]") {
  // A run of Ns reverse-complements to a run of Ns of the same length. Edge
  // case for ambiguous-IUPAC-bearing sequences from low-confidence regions.
  CHECK(RevComp(std::string_view("NNNNN")) == "NNNNN");
  CHECK(RevComp(std::string_view("nnn")) == "nnn");
}

TEST_CASE("RevComp(string_view) preserves mixed case", "[lancet][base][RevComp]") {
  // Mixed-case input from soft-masked aligners; the implementation must not
  // upcase or downcase as a side effect.
  CHECK(RevComp(std::string_view("AcGt")) == "aCgT");
  CHECK(RevComp(std::string_view("aTGc")) == "gCAt");
}

TEST_CASE("RevComp(string_view) handles a fully lowercase sequence", "[lancet][base][RevComp]") {
  CHECK(RevComp(std::string_view("acgt")) == "acgt");
  CHECK(RevComp(std::string_view("gattaca")) == "tgtaatc");
}

TEST_CASE("RevComp is involutive on canonical DNA (RevComp(RevComp(x)) == x)",
          "[lancet][base][RevComp]") {
  // The reverse-complement is its own inverse (involution): applying it
  // twice gets you back to the original sequence. A regression that broke
  // either the reversal or the complement (e.g. complement-then-reverse vs
  // reverse-then-complement bug) would surface here.
  for (auto const& seq :
       {std::string_view("A"), std::string_view("ACGT"), std::string_view("GATTACA"),
        std::string_view("AcGtNn"), std::string_view(""), std::string_view("NNNNN")}) {
    INFO("seq=\"" << seq << "\"");
    CHECK(RevComp(RevComp(seq)) == seq);
  }
}

}  // namespace lancet::base::tests
