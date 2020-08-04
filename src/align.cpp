#include "lancet/align.h"

#include <algorithm>
#include <cstdint>
#include <utility>

#include "absl/strings/str_format.h"

namespace lancet {
// Original:
static constexpr int MATCH_SCORE = 2;
static constexpr int MISMATCH_SCORE = -4;
static constexpr int GAP_OPEN_SCORE = -8;
static constexpr int GAP_EXTEND_SCORE = -1;

inline auto Compare(const char& s, const char& t) -> int { return s == t ? MATCH_SCORE : MISMATCH_SCORE; }

struct Cell {
  explicit Cell(int s = 0, char tb = '*') : score(s), trace_back(tb) {}

  int score;        // NOLINT
  char trace_back;  // NOLINT
};

static inline auto MaxScoreX(int a, int b) -> Cell { return a > b ? Cell(a, '-') : Cell(b, '<'); }
static inline auto MaxScoreY(int a, int b) -> Cell { return a > b ? Cell(a, '|') : Cell(b, '^'); }
static inline auto MaxScoreXY(const Cell& scorex, const Cell& scorey, int scorez) -> Cell {
  auto result = Cell(scorez, '\\');

  if (scorex.score > result.score) {
    result.score = scorex.score;
    result.trace_back = '<';
  }

  if (scorey.score > result.score) {
    result.score = scorey.score;
    result.trace_back = '^';
  }

  return result;
}

using AlignedBases = std::pair<char, char>;
auto Align(std::string_view ref, std::string_view qry) -> AlignedSequences {
  const auto refLen = ref.length();
  const auto qryLen = qry.length();

  auto x = detail::NaiveTwoDimArray<Cell>{refLen + 2, qryLen + 2};
  auto y = detail::NaiveTwoDimArray<Cell>{refLen + 2, qryLen + 2};
  auto m = detail::NaiveTwoDimArray<Cell>{refLen + 2, qryLen + 2};

  for (std::size_t j = 0; j <= qryLen; ++j) {
    x.at(0, j).score = GAP_OPEN_SCORE + (static_cast<int>(j) * GAP_EXTEND_SCORE);
    x.at(0, j).trace_back = '^';
    m.at(0, j) = x.at(0, j);
  }

  for (std::size_t i = 0; i <= refLen; ++i) {
    y.at(i, 0).score = GAP_OPEN_SCORE + (static_cast<int>(i) * GAP_EXTEND_SCORE);
    y.at(i, 0).trace_back = '<';
    m.at(i, 0) = y.at(i, 0);
  }

  // qi -> query sequence index
  // ri -> reference sequence index
  for (std::size_t qi = 1; qi <= qryLen; qi++) {
    for (std::size_t ri = 1; ri <= refLen; ri++) {
      x.at(ri, qi) = MaxScoreX(x.at(ri - 1, qi).score + GAP_EXTEND_SCORE, m.at(ri - 1, qi).score + GAP_OPEN_SCORE);
      y.at(ri, qi) = MaxScoreY(y.at(ri, qi - 1).score + GAP_EXTEND_SCORE, m.at(ri, qi - 1).score + GAP_OPEN_SCORE);
      m.at(ri, qi) =
          MaxScoreXY(x.at(ri, qi), y.at(ri, qi), m.at(ri - 1, qi - 1).score + Compare(ref[ri - 1], qry[qi - 1]));
    }
  }

  AlignedBases alnBases;
  std::int64_t i = refLen;
  std::int64_t j = qryLen;
  bool forcey = false;
  bool forcex = false;

  std::string refAln;
  std::string qryAln;
  refAln.reserve(std::max(refLen, qryLen) + std::min(refLen, qryLen));
  qryAln.reserve(std::max(refLen, qryLen) + std::min(refLen, qryLen));

  /*
   * gdb --args ~/dev/projects/Lancet/lancet/build-gcc9-debug/bin/lancet --tumor
   * ~/dev/projects/Lancet/VirtualTumors/Tumor_INDEL.NA12892_mem_binomial_indel.final.bam --normal
   * ~/dev/projects/Lancet/VirtualTumors/Normal.NA12892_mem_normal.final.bam --reference
   * ~/dev/projects/Lancet/human_g1k_v37_decoy.fasta --out-vcf indel.r1.vcf.gz --log-level TRACE --region
   * 1:3553050-3553150
   *
   * catch signal
   * catch catch
   *
   *
   * gdb --args /nfs/sw/lancet/lancet-1.0.7/lancet --tumor
   * ~/dev/projects/Lancet/VirtualTumors/Tumor_INDEL.NA12892_mem_binomial_indel.final.bam --normal
   * ~/dev/projects/Lancet/VirtualTumors/Normal.NA12892_mem_normal.final.bam --ref
   * ~/dev/projects/Lancet/human_g1k_v37_decoy.fasta --reg 1:3553050-3553150 --verbose
   *
   */

  while (i > 0 || j > 0) {
    char t = m.at(i, j).trace_back;
    if (t == '*') break;

    if (forcex) {
      alnBases.first = ref[i - 1];
      alnBases.second = ALIGNMENT_GAP;
      if (x.at(i, j).trace_back == '<') forcex = false;
      --i;
    } else if (t == '<') {
      alnBases.first = ref[i - 1];
      alnBases.second = ALIGNMENT_GAP;
      if (x.at(i, j).trace_back == '-') forcex = true;
      --i;
    } else if (forcey) {
      alnBases.first = ALIGNMENT_GAP;
      alnBases.second = qry[j - 1];
      if (y.at(i, j).trace_back == '^') forcey = false;
      --j;
    } else if (t == '^') {
      alnBases.first = ALIGNMENT_GAP;
      alnBases.second = qry[j - 1];
      if (y.at(i, j).trace_back == '|') forcey = true;
      --j;
    } else if (t == '\\') {
      alnBases.first = ref[i - 1];
      alnBases.second = qry[j - 1];
      --i;
      --j;
    } else {
      throw std::runtime_error{
          absl::StrFormat("unexpected error in SW alignment. ref: %s, query: %s", ref, qry)};
    }

    refAln.push_back(alnBases.first);
    qryAln.push_back(alnBases.second);
  }

  std::reverse(refAln.begin(), refAln.end());
  std::reverse(qryAln.begin(), qryAln.end());
  return AlignedSequences{refAln, qryAln};
}
}  // namespace lancet
