#include "lancet/hts/cigar_utils.h"

#include <vector>

#include "catch_amalgamated.hpp"

extern "C" {
#include "htslib/sam.h"
}

#include "lancet/base/types.h"
#include "lancet/hts/cigar_unit.h"

using namespace lancet::hts;

namespace {

/// Helper: build a CigarUnit from length + BAM_CIGAR op code (e.g. BAM_CMATCH).
[[nodiscard]] auto MakeCigar(u32 len, u32 bam_op) -> CigarUnit {
  return CigarUnit(bam_cigar_gen(len, bam_op));
}

/// Helper: encode a DNA string into numeric 0-4 values matching the genotyper
/// encoding convention (A=0, C=1, G=2, T=3, N/other=4).
[[nodiscard]] auto EncodeSeq(std::string_view seq) -> std::vector<u8> {
  // Same ENCODE_TABLE as genotyper_detail but inline for test isolation.
  // NOLINTNEXTLINE(readability-magic-numbers)
  static constexpr std::array<u8, 256> ENCODE = []() {
    std::array<u8, 256> table{};
    table.fill(4);
    table['A'] = 0; table['a'] = 0;
    table['C'] = 1; table['c'] = 1;
    table['G'] = 2; table['g'] = 2;
    table['T'] = 3; table['t'] = 3;
    return table;
  }();

  std::vector<u8> result(seq.size());
  for (usize i = 0; i < seq.size(); ++i) {
    result[i] = ENCODE[static_cast<u8>(seq[i])];
  }
  return result;
}

}  // namespace

// ============================================================================
// ComputeEditDistance tests
// ============================================================================
TEST_CASE("ComputeEditDistance: perfect match has NM=0", "[lancet][hts][CigarUtils]") {
  // CIGAR: 10M, query == target → NM = 0
  const std::vector<CigarUnit> cigar{MakeCigar(10, BAM_CMATCH)};
  const auto query = EncodeSeq("ACGTACGTAC");
  const auto target = EncodeSeq("ACGTACGTAC");
  CHECK(ComputeEditDistance(cigar, query, target) == 0);
}

TEST_CASE("ComputeEditDistance: mismatches under M ops", "[lancet][hts][CigarUtils]") {
  // CIGAR: 5M, 2 mismatches at positions 1 and 3
  const std::vector<CigarUnit> cigar{MakeCigar(5, BAM_CMATCH)};
  const auto query  = EncodeSeq("ATGTA");
  const auto target = EncodeSeq("AAGAA");
  // Differences: pos1 (T vs A), pos3 (T vs A) → NM = 2
  CHECK(ComputeEditDistance(cigar, query, target) == 2);
}

TEST_CASE("ComputeEditDistance: insertion adds to NM", "[lancet][hts][CigarUtils]") {
  // CIGAR: 3M 2I 3M → query has 8 bases, target has 6 bases
  // NM includes 2 insertion bases
  const std::vector<CigarUnit> cigar{
      MakeCigar(3, BAM_CMATCH), MakeCigar(2, BAM_CINS), MakeCigar(3, BAM_CMATCH)};
  const auto query  = EncodeSeq("ACGTTACG");  // 8 bases
  const auto target = EncodeSeq("ACGACG");     // 6 bases
  // M:ACG matches ACG (0), I:TT (+2), M:ACG matches ACG (0) → NM = 2
  CHECK(ComputeEditDistance(cigar, query, target) == 2);
}

TEST_CASE("ComputeEditDistance: deletion adds to NM", "[lancet][hts][CigarUtils]") {
  // CIGAR: 3M 2D 3M → query has 6 bases, target has 8 bases
  const std::vector<CigarUnit> cigar{
      MakeCigar(3, BAM_CMATCH), MakeCigar(2, BAM_CDEL), MakeCigar(3, BAM_CMATCH)};
  const auto query  = EncodeSeq("ACGACG");     // 6 bases
  const auto target = EncodeSeq("ACGTTACG");   // 8 bases
  // M:ACG (0), D:TT (+2), M:ACG (0) → NM = 2
  CHECK(ComputeEditDistance(cigar, query, target) == 2);
}

TEST_CASE("ComputeEditDistance: soft clips excluded from NM", "[lancet][hts][CigarUtils]") {
  // CIGAR: 3S 5M 2S → only the 5M portion counts
  const std::vector<CigarUnit> cigar{
      MakeCigar(3, BAM_CSOFT_CLIP), MakeCigar(5, BAM_CMATCH), MakeCigar(2, BAM_CSOFT_CLIP)};
  const auto query  = EncodeSeq("NNNACGTANN");  // 10 bases (3S + 5M + 2S)
  const auto target = EncodeSeq("ACGTA");       // 5 bases (only matched region)
  CHECK(ComputeEditDistance(cigar, query, target) == 0);
}

TEST_CASE("ComputeEditDistance: mixed CIGAR with mismatches, ins, del", "[lancet][hts][CigarUtils]") {
  // CIGAR: 4M 1I 2M 1D 4M
  // Ref:   A C G T - A C G T A A A  (11 consumed from target)
  // Query: A C A T C A C - T A A A  (11 consumed from query)
  // Mismatches under M: pos2 (A vs G) = 1 mismatch
  // I: 1 base, D: 1 base → NM = 1 + 1 + 1 = 3
  const std::vector<CigarUnit> cigar{
      MakeCigar(4, BAM_CMATCH), MakeCigar(1, BAM_CINS),
      MakeCigar(2, BAM_CMATCH), MakeCigar(1, BAM_CDEL),
      MakeCigar(4, BAM_CMATCH)};
  const auto query  = EncodeSeq("ACATCACTAAA");  // 11 bases
  const auto target = EncodeSeq("ACGTACGTAAA");  // 11 bases
  CHECK(ComputeEditDistance(cigar, query, target) == 3);
}

TEST_CASE("ComputeEditDistance: =/X operators work correctly", "[lancet][hts][CigarUtils]") {
  // CIGAR: 3= 1X 2= → 1 mismatch from X, 0 from =
  const std::vector<CigarUnit> cigar{
      MakeCigar(3, BAM_CEQUAL), MakeCigar(1, BAM_CDIFF), MakeCigar(2, BAM_CEQUAL)};
  const auto query  = EncodeSeq("ACGTAC");
  const auto target = EncodeSeq("ACGAAC");
  CHECK(ComputeEditDistance(cigar, query, target) == 1);
}

// ============================================================================
// CigarRefPosToQueryPos tests
// ============================================================================
TEST_CASE("CigarRefPosToQueryPos: simple match", "[lancet][hts][CigarUtils]") {
  // 10M: ref pos 5 → query pos 5
  const std::vector<CigarUnit> cigar{MakeCigar(10, BAM_CMATCH)};
  CHECK(CigarRefPosToQueryPos(cigar, 0) == 0);
  CHECK(CigarRefPosToQueryPos(cigar, 5) == 5);
  CHECK(CigarRefPosToQueryPos(cigar, 9) == 9);
}

TEST_CASE("CigarRefPosToQueryPos: insertion shifts query ahead", "[lancet][hts][CigarUtils]") {
  // 3M 2I 5M: after the insertion, query is 2 bases ahead of ref
  // Ref pos 0-2 → query 0-2 (through first 3M)
  // Then 2I: query advances to 5, ref stays at 3
  // Ref pos 3 → query pos 5
  const std::vector<CigarUnit> cigar{
      MakeCigar(3, BAM_CMATCH), MakeCigar(2, BAM_CINS), MakeCigar(5, BAM_CMATCH)};
  CHECK(CigarRefPosToQueryPos(cigar, 2) == 2);
  CHECK(CigarRefPosToQueryPos(cigar, 3) == 5);
  CHECK(CigarRefPosToQueryPos(cigar, 7) == 9);
}

TEST_CASE("CigarRefPosToQueryPos: deletion maps to query at del start", "[lancet][hts][CigarUtils]") {
  // 3M 2D 5M: ref pos 3 falls inside the deletion
  // No query base at ref pos 3, so return qpos at deletion start (= 3)
  const std::vector<CigarUnit> cigar{
      MakeCigar(3, BAM_CMATCH), MakeCigar(2, BAM_CDEL), MakeCigar(5, BAM_CMATCH)};
  CHECK(CigarRefPosToQueryPos(cigar, 2) == 2);
  CHECK(CigarRefPosToQueryPos(cigar, 3) == 3);  // inside deletion
  CHECK(CigarRefPosToQueryPos(cigar, 4) == 3);  // still inside deletion
  CHECK(CigarRefPosToQueryPos(cigar, 5) == 3);  // first base after deletion
  CHECK(CigarRefPosToQueryPos(cigar, 6) == 4);
}

TEST_CASE("CigarRefPosToQueryPos: soft clip advances query before match", "[lancet][hts][CigarUtils]") {
  // 3S 5M: soft clip doesn't consume reference but advances query
  // Ref pos 0 corresponds to query pos 3
  const std::vector<CigarUnit> cigar{
      MakeCigar(3, BAM_CSOFT_CLIP), MakeCigar(5, BAM_CMATCH)};
  CHECK(CigarRefPosToQueryPos(cigar, 0) == 3);
  CHECK(CigarRefPosToQueryPos(cigar, 4) == 7);
}
