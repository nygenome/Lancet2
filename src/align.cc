#include "lancet2/align.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "absl/strings/str_format.h"

namespace {
template <typename T>
class TwoDimVector {
 public:
  TwoDimVector(std::size_t nrows, std::size_t ncols) : numRows(nrows), numCols(ncols), data(nrows * ncols) {}
  TwoDimVector() = delete;

  [[nodiscard]] auto at(std::size_t row, std::size_t col) -> T& {
    check(row, col);
    return data.at(row * numCols + col);
  }

  [[nodiscard]] auto at(std::size_t row, std::size_t col) const -> const T& {
    check(row, col);
    return data.at(row * numCols + col);
  }

 private:
  std::size_t numRows;
  std::size_t numCols;
  std::vector<T> data;  // NOLINT

  auto check(std::size_t row, std::size_t col) -> void {
    if ((row >= numRows) or (col >= numCols)) {
      throw std::out_of_range(
          absl::StrFormat("requestedRow=%d, numRows=%d, requestedCol=%d, numCols=%d", row, numRows, col, numCols));
    }
  }
};
}  // namespace

namespace lancet2 {
// Original:
static constexpr int MATCH_SCORE = 2;
static constexpr int MISMATCH_SCORE = -4;
static constexpr int GAP_OPEN_SCORE = -8;
static constexpr int GAP_EXTEND_SCORE = -1;

inline auto Compare(const char& s, const char& t) -> int { return s == t ? MATCH_SCORE : MISMATCH_SCORE; }

struct Cell {
  Cell() = default;
  Cell(int s, char tb) : score(s), trace_back(tb) {}

  int score = 0;          // NOLINT
  char trace_back = '*';  // NOLINT
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

static constexpr inline auto IsValidBase(const char& b) -> bool { return b == 'A' || b == 'C' || b == 'G' || b == 'T'; }

using AlignedBases = std::pair<char, char>;
auto Align(std::string_view ref, std::string_view qry) -> AlignedSequences {
  const auto refLen = ref.length();
  const auto qryLen = qry.length();

  TwoDimVector<Cell> x(refLen + 2, qryLen + 2);
  TwoDimVector<Cell> y(refLen + 2, qryLen + 2);
  TwoDimVector<Cell> m(refLen + 2, qryLen + 2);

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
      const auto cmpScore = Compare(ref[ri - 1], qry[qi - 1]);
      x.at(ri, qi) = MaxScoreX(x.at(ri - 1, qi).score + GAP_EXTEND_SCORE, m.at(ri - 1, qi).score + GAP_OPEN_SCORE);
      y.at(ri, qi) = MaxScoreY(y.at(ri, qi - 1).score + GAP_EXTEND_SCORE, m.at(ri, qi - 1).score + GAP_OPEN_SCORE);
      m.at(ri, qi) = MaxScoreXY(x.at(ri, qi), y.at(ri, qi), m.at(ri - 1, qi - 1).score + cmpScore);
    }
  }

  AlignedBases alnBases;
  int i = static_cast<int>(refLen);
  int j = static_cast<int>(qryLen);
  bool forcey = false;
  bool forcex = false;

  std::string refAln;
  std::string qryAln;
  refAln.reserve(std::max(refLen, qryLen) + std::min(refLen, qryLen));
  qryAln.reserve(std::max(refLen, qryLen) + std::min(refLen, qryLen));

  while (i > 0 || j > 0) {
    char t = m.at(i, j).trace_back;
    if (t == '*') break;

    if (forcex) {
      alnBases.first = ref.at(i - 1);
      alnBases.second = ALIGN_GAP;
      if (x.at(i, j).trace_back == '<') forcex = false;
      --i;
    } else if (t == '<') {
      alnBases.first = ref.at(i - 1);
      alnBases.second = ALIGN_GAP;
      if (x.at(i, j).trace_back == '-') forcex = true;
      --i;
    } else if (forcey) {
      alnBases.first = ALIGN_GAP;
      alnBases.second = qry.at(j - 1);
      if (y.at(i, j).trace_back == '^') forcey = false;
      --j;
    } else if (t == '^') {
      alnBases.first = ALIGN_GAP;
      alnBases.second = qry.at(j - 1);
      if (y.at(i, j).trace_back == '|') forcey = true;
      --j;
    } else if (t == '\\') {
      alnBases.first = ref.at(i - 1);
      alnBases.second = qry.at(j - 1);
      --i;
      --j;
    } else {
      throw std::runtime_error{absl::StrFormat("unexpected error in SW alignment. ref: %s, query: %s", ref, qry)};
    }

    refAln.push_back(IsValidBase(alnBases.first) ? alnBases.first : ALIGN_GAP);
    qryAln.push_back(IsValidBase(alnBases.second) ? alnBases.second : ALIGN_GAP);
  }

  std::reverse(refAln.begin(), refAln.end());
  std::reverse(qryAln.begin(), qryAln.end());
  return AlignedSequences{refAln, qryAln};
}

auto TrimEndGaps(AlignedSequencesView* aln) -> std::size_t {
  // Trim end GAPS and adjust end alignments until both ends in ref and qry have no GAPS
  std::size_t refStartTrim = 0;
  std::size_t start = 0;
  std::size_t end = aln->ref.length() - 1;

  const auto startGap = aln->ref[start] == ALIGN_GAP || aln->qry[start] == ALIGN_GAP;
  const auto endGap = aln->ref[end] == ALIGN_GAP || aln->qry[end] == ALIGN_GAP;

  if (startGap || endGap) {
    // move start until no begin alignment gaps are found
    while (aln->ref[start] == ALIGN_GAP || aln->qry[start] == ALIGN_GAP) {
      if (aln->ref[start] == ALIGN_GAP) refStartTrim++;
      start++;
    }

    // move end until no end alignment gaps are found
    while (aln->ref[end] == ALIGN_GAP || aln->qry[end] == ALIGN_GAP) {
      end--;
    }

    aln->ref = aln->ref.substr(start, end - start);
    aln->qry = aln->qry.substr(start, end - start);
  }

  return refStartTrim;
}
}  // namespace lancet2
