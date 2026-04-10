#include "lancet/base/longdust_scorer.h"

#include "lancet/base/types.h"
#include "lancet/hts/reference.h"

#include "catch_amalgamated.hpp"
#include "lancet_test_config.h"

#include <array>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <string_view>

#include <cmath>

// longdust.c and kalloc.c are compiled into TestLancet2 via FetchContent.
// This header provides access to longdust's internal f[] table and Q-score
// functions for exact cross-validation against our C++ LongdustQScorer.
#include "base/longdust_test_helpers.h"

namespace {

// ── Shared constants & helpers ──────────────────────────────────────────────

constexpr auto CHM13_REF_NAME = "chm13v2.0.fa.gz";
constexpr auto CALIBRATION_TSV = "longdust_calibration_test_regions.tsv";

/// Build a tandem repeat string: motif repeated `copies` times.
inline auto BuildRepeat(std::string_view motif, usize copies) -> std::string {
  std::string result;
  result.reserve(motif.size() * copies);
  for (usize iter = 0; iter < copies; ++iter) {
    result.append(motif);
  }
  return result;
}

/// Generate a pseudorandom DNA sequence of given length.
inline auto RandomDna(usize length, u64 seed = 42) -> std::string {
  static constexpr std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<usize> dist(0, 3);
  std::string result;
  result.reserve(length);
  for (usize iter = 0; iter < length; ++iter) {
    result.push_back(BASES.at(dist(gen)));
  }
  return result;
}

/// RAII wrapper for longdust data (auto-destroys on scope exit).
struct LdContext {
  ld_opt_t mOpt{};
  ld_data_t* mData = nullptr;

  // NOLINTBEGIN(cppcoreguidelines-prefer-member-initializer)
  LdContext() {
    ld_opt_init(&mOpt);
    mData = ld_data_init(nullptr, &mOpt);
  }
  // NOLINTEND(cppcoreguidelines-prefer-member-initializer)
  ~LdContext() { ld_data_destroy(mData); }
  LdContext(LdContext const&) = delete;
  auto operator=(LdContext const&) -> LdContext& = delete;
  LdContext(LdContext&&) = delete;
  auto operator=(LdContext&&) -> LdContext& = delete;
};

/// Row from the calibration TSV (new format with BED coordinates and pre-computed scores).
struct CalibrationRow {
  std::string mChrom;   // e.g. "chr1"
  i64 mStart = 0;       // 0-based start
  i64 mEnd = 0;         // 0-based end (exclusive)
  std::string mSource;  // category (e.g. "GIAB/LowComplexity")
  std::string mName;    // subcategory (e.g. "SimpleRepeat_diTR_11to50_slop5")
  std::string mRegion;  // samtools region (e.g. "chr1:100-200")
  int mScale = 0;
  int mRegionLength = 0;
  double mExpectedScore{};  // pre-computed LCR score
};

/// Parse the calibration TSV into a vector of rows.
inline auto LoadCalibrationTsv(std::filesystem::path const& path) -> std::vector<CalibrationRow> {
  std::ifstream infile(path);
  std::vector<CalibrationRow> rows;
  std::string line;
  std::getline(infile, line);  // skip header
  while (std::getline(infile, line)) {
    if (line.empty()) continue;
    CalibrationRow row;
    std::istringstream iss(line);
    std::string start_s;
    std::string end_s;
    std::string scale_s;
    std::string rlen_s;
    std::string score_s;
    std::getline(iss, row.mChrom, '\t');
    std::getline(iss, start_s, '\t');
    std::getline(iss, end_s, '\t');
    std::getline(iss, row.mSource, '\t');
    std::getline(iss, row.mName, '\t');
    std::getline(iss, row.mRegion, '\t');
    std::getline(iss, scale_s, '\t');
    std::getline(iss, rlen_s, '\t');
    std::getline(iss, score_s, '\t');
    row.mStart = std::stol(start_s);
    row.mEnd = std::stol(end_s);
    row.mScale = std::stoi(scale_s);
    row.mRegionLength = std::stoi(rlen_s);
    row.mExpectedScore = score_s.empty() ? 0.0 : std::stod(score_s);
    rows.push_back(std::move(row));
  }
  return rows;
}

}  // namespace

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PART 1: Longdust Cross-Validation                                       ║
// ║  Verify exact mathematical parity between LongdustQScorer and longdust C.║
// ╚══════════════════════════════════════════════════════════════════════════╝

// ============================================================================
// 1a. Real CHM13 sequences — score parity on authentic genomic DNA
// ============================================================================

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cross-validation: LongdustQScorer vs longdust on real CHM13 sequences",
          "[lancet][base][longdust_scorer][cross_validation]") {
  auto const ref_path = MakePath(FULL_DATA_DIR, CHM13_REF_NAME);
  REQUIRE(std::filesystem::exists(ref_path));
  lancet::hts::Reference const ref(ref_path);
  lancet::base::LongdustQScorer const scorer(7, 5001, 0.5);  // gc_frac=0.5 for longdust parity
  LdContext lctx;
  REQUIRE(lctx.mData != nullptr);
  constexpr double EPS = 1e-9;

  // Telomere: (TTAGGG)n — highest complexity score
  // Source: chm13v2.0_telomere.bed chr1:0-3000
  SECTION("chr1 p-arm telomere (200bp)") {
    auto const region = ref.MakeRegion("chr1:1-200");
    auto const seq = std::string(region.SeqView());
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 200)).epsilon(EPS));
    CHECK(scorer.Score(seq) > 0.6);
  }

  // Source: chm13v2.0_telomere.bed chr2:242693800-242696752
  SECTION("chr2 q-arm telomere (200bp)") {
    auto const region = ref.MakeRegion("chr2:242696553-242696752");
    auto const seq = std::string(region.SeqView());
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 200)).epsilon(EPS));
    CHECK(scorer.Score(seq) > 0.5);
  }

  // HSAT2: ~23bp pericentromeric satellite
  // Source: chm13v2.0_censat_v2.1.bed chr1:126828704-126838321 hsat2_1_1(B)
  SECTION("chr1 HSAT2 pericentromeric (200bp)") {
    auto const region = ref.MakeRegion("chr1:126828705-126828904");
    auto const seq = std::string(region.SeqView());
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 200)).epsilon(EPS));
    CHECK(scorer.Score(seq) > 0.3);
  }

  // Alpha satellite HOR: ~171bp centromeric monomer
  // Source: chm13v2.0_censat_v2.1.bed chr1:121619169-121625213 hor_1_1
  SECTION("chr1 alpha-sat HOR (300bp)") {
    auto const region = ref.MakeRegion("chr1:121619170-121619469");
    auto const seq = std::string(region.SeqView());
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 300)).epsilon(EPS));
    CHECK(scorer.Score(seq) > 0.0);
  }

  // Ajax (AATGG)n pentanucleotide satellite
  // Source: chm13v2.0_new-satellites_2022DEC.bed chr1:153734-153878 ajax
  SECTION("chr1 ajax satellite (144bp)") {
    auto const region = ref.MakeRegion("chr1:153735-153878");
    auto const seq = std::string(region.SeqView());
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 144)).epsilon(EPS));
    CHECK(scorer.Score(seq) > 0.2);
  }

  // Non-repetitive unique region (control)
  SECTION("chr1 non-repetitive unique (200bp)") {
    auto const region = ref.MakeRegion("chr1:10000001-10000200");
    auto const seq = std::string(region.SeqView());
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 200)).epsilon(EPS));
    CHECK(scorer.Score(seq) < 0.3);
  }

  // Subtelomeric region with interspersed simple repeats
  SECTION("chr1 subtelomeric (200bp)") {
    auto const region = ref.MakeRegion("chr1:3001-3200");
    auto const seq = std::string(region.SeqView());
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 200)).epsilon(EPS));
  }
}

// ============================================================================
// 1b. f(ℓ) table parity — internal Poisson expectation function
// ============================================================================

TEST_CASE("Cross-validation: f(ℓ) table matches longdust",
          "[lancet][base][longdust_scorer][cross_validation]") {
  LdContext const lctx;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* ldat = reinterpret_cast<ld_test_data_s const*>(lctx.mData);
  lancet::base::LongdustQScorer const scorer(7, 5001, 0.5);  // gc_frac=0.5 for longdust parity

  // For a homopolymer of ℓ+k-1 bases: Q = lgamma(ℓ+1) - f(ℓ)
  // so f(ℓ) = lgamma(ℓ+1) - ScoreOneStrand() * ℓ
  for (int const ell : {10, 50, 100, 200, 500, 1000, 2000, 4000}) {
    auto const seq = std::string(ell + 6, 'A');  // k=7 → ℓ k-mers
    double const ld_f = ldat->f[ell];
    double const our_f =
        std::lgamma(static_cast<f64>(ell + 1)) - (scorer.ScoreOneStrand(seq) * ell);
    INFO("ℓ=" << ell << " ld_f=" << ld_f << " our_f=" << our_f);
    CHECK(our_f == Catch::Approx(ld_f).epsilon(1e-9));
  }
}

// ============================================================================
// 1c. Synthetic repeats — exact tandem repeat formula validation
// ============================================================================

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cross-validation: LongdustQScorer vs longdust on synthetic repeats",
          "[lancet][base][longdust_scorer][cross_validation]") {
  lancet::base::LongdustQScorer const scorer(7, 5001, 0.5);  // gc_frac=0.5 for longdust parity
  LdContext const lctx;
  constexpr double EPS = 1e-9;

  SECTION("Homopolymer poly-A 50bp") {
    auto const seq = std::string(50, 'A');
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 50)).epsilon(EPS));
  }

  SECTION("Homopolymer poly-A 100bp") {
    auto const seq = std::string(100, 'A');
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 100)).epsilon(EPS));
  }

  SECTION("Dinucleotide (CA)x50") {
    auto const seq = BuildRepeat("CA", 50);
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 100)).epsilon(EPS));
  }

  SECTION("Trinucleotide (CAG)x40") {
    auto const seq = BuildRepeat("CAG", 40);
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 120)).epsilon(EPS));
  }

  SECTION("Random DNA 200bp") {
    auto const seq = RandomDna(200);
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 200)).epsilon(EPS));
  }

  SECTION("Random DNA 500bp") {
    auto const seq = RandomDna(500, 789);
    CHECK(scorer.Score(seq) ==
          Catch::Approx(ld_test_q_both_strands(lctx.mData, seq.c_str(), 500)).epsilon(EPS));
  }
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PART 2: Scorer Properties                                               ║
// ║  Unit tests for edge cases, monotonicity, and strand behaviour.          ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("Score: zero for short, empty, or N-only sequences", "[lancet][base][longdust_scorer]") {
  lancet::base::LongdustQScorer const scorer(7);
  CHECK(scorer.Score("") == 0.0);
  CHECK(scorer.Score("ATCG") == 0.0);
  CHECK(scorer.Score("ATCGAT") == 0.0);
  CHECK(scorer.Score("NNNNNNNNNNNNNNNNNNN") == 0.0);
}

TEST_CASE("Score: near-zero for random DNA", "[lancet][base][longdust_scorer]") {
  lancet::base::LongdustQScorer const scorer(7);
  CHECK(scorer.Score(RandomDna(100)) < 0.1);
  CHECK(scorer.Score(RandomDna(200, 123)) < 0.1);
  CHECK(scorer.Score(RandomDna(500, 456)) < 0.1);
}

TEST_CASE("Score: detects homopolymer runs", "[lancet][base][longdust_scorer]") {
  lancet::base::LongdustQScorer const scorer(7);
  CHECK(scorer.Score(std::string(10, 'A')) > 0.0);
  CHECK(scorer.Score(std::string(20, 'A')) > 0.6);
  CHECK(scorer.Score(std::string(50, 'A')) > 1.0);
  CHECK(scorer.Score(std::string(20, 'A')) < scorer.Score(std::string(50, 'A')));
  CHECK(scorer.Score(std::string(50, 'A')) < scorer.Score(std::string(100, 'A')));
}

TEST_CASE("Score: increases with repeat copy number", "[lancet][base][longdust_scorer]") {
  lancet::base::LongdustQScorer const scorer(7);
  auto const score5 = scorer.Score(BuildRepeat("TTAGGG", 5));
  auto const score10 = scorer.Score(BuildRepeat("TTAGGG", 10));
  auto const score20 = scorer.Score(BuildRepeat("TTAGGG", 20));
  CHECK(score5 < score10);
  CHECK(score10 <= score20);
}

TEST_CASE("Score: uses both strands (max of fwd/rev)", "[lancet][base][longdust_scorer]") {
  lancet::base::LongdustQScorer const scorer(7);
  auto const poly_t = std::string(30, 'T');
  CHECK(scorer.Score(poly_t) > 0.6);
  CHECK(scorer.Score(poly_t) >= scorer.ScoreOneStrand(poly_t));
}

TEST_CASE("Score: case-insensitive", "[lancet][base][longdust_scorer]") {
  lancet::base::LongdustQScorer const scorer(7);
  auto const upper = BuildRepeat("TTAGGG", 20);
  std::string lower = upper;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  CHECK(scorer.Score(upper) == Catch::Approx(scorer.Score(lower)));
}

TEST_CASE("Score: N bases reduce score", "[lancet][base][longdust_scorer]") {
  lancet::base::LongdustQScorer const scorer(7);
  CHECK(scorer.Score("AAAAAAANAAAAAAA") < scorer.Score(std::string(15, 'A')));
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PART 3: Score Formatting & Constants                                    ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("FormatComplexityScore: precision and trailing-zero stripping",
          "[lancet][base][longdust_scorer]") {
  using lancet::base::FormatComplexityScore;
  CHECK(FormatComplexityScore(0.0) == "0");
  CHECK(FormatComplexityScore(1.0) == "1");
  CHECK(FormatComplexityScore(0.5) == "0.5");
  CHECK(FormatComplexityScore(0.123) == "0.123");
  CHECK(FormatComplexityScore(0.1236) == "0.124");
  CHECK(FormatComplexityScore(1.9999) == "2");
  CHECK(FormatComplexityScore(0.100) == "0.1");
}

TEST_CASE("FormatComplexityScores: comma-separated output", "[lancet][base][longdust_scorer]") {
  using lancet::base::FormatComplexityScores;
  CHECK(FormatComplexityScores(std::initializer_list<f64>{0, 0, 0, 0, 0}) == "0,0,0,0,0");
  CHECK(FormatComplexityScores(std::initializer_list<f64>{1.8, 0.7, 0.3, 0.15, 0.1}) ==
        "1.8,0.7,0.3,0.15,0.1");
}

TEST_CASE("Longdust k-mer size constants", "[lancet][base][longdust_scorer]") {
  CHECK(lancet::base::LONGDUST_FLANK_K == 4);
  CHECK(lancet::base::LONGDUST_HAPLOTYPE_K == 7);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PART 4: Multi-Scale Calibration from TSV                                ║
// ║                                                                          ║
// ║  Validates pre-computed LCR scores from the test calibration TSV.        ║
// ║  The TSV is generated by ScoreRegionLCR with balanced sampling           ║
// ║  (≤5000 regions per subcategory, ~8 sampled for test TSV).               ║
// ║                                                                          ║
// ║  TSV format: chrom start end source name region scale region_length score║
// ║  Generated by: ScoreRegionLCR (tests/base/score_region_lcr.cpp)          ║
// ╚══════════════════════════════════════════════════════════════════════════╝
//
// CALIBRATION SUMMARY — from CHM13v2.0 balanced analysis (133,714 regions)
// ─────────────────────────────────────────────────────────────────────────
// Balanced: ≤5000 per subcategory, 135 subcategories, 4 scales
// Full data: data/lcr_calibration_balanced_scored.tsv.gz (534,856 rows)
//            data/lcr_analysis_by_category.txt (category breakdown)
//            data/lcr_analysis_by_subcategory.txt (subcategory breakdown)
//
//   CATEGORY-LEVEL (median scores):
//
//   Source Category       100bp   200bp  1000bp  2000bp
//   ─────────────────── ─────── ─────── ─────── ───────
//   Telomere (TTAGGG)n    1.897   2.554   3.678   4.037
//   NewSatellite          0.138   0.286   0.273   0.498
//   GIAB/LowComplexity    0.236   0.138   0.081   0.094
//   CenSat                0.005   0.017   0.113   0.174
//   RepeatMasker (all)    0.005   0.017   0.052   0.079
//   GIAB/Union            0.000   0.010   0.048   0.073
//   SegDups               0.005   0.010   0.039   0.065
//
//   GIAB SUBCATEGORIES (median scores, sorted by 2000bp desc):
//   triTR_gt200             2.527   3.212   0.929   0.529
//   quadTR_gt200            2.262   2.937   0.755   0.386
//   satellites              0.013   0.037   0.240   0.375
//   triTR_51to200           1.295   0.785   0.288   0.198
//   quadTR_51to200          1.198   0.804   0.245   0.179
//   diTR_51to200            1.198   0.595   0.155   0.121
//   homopolymer_gt20        0.399   0.201   0.088   0.094
//   homopolymer_4to6        0.000   0.010   0.038   0.064
//   SegDups                 0.005   0.010   0.039   0.065
//
//   REPEATMASKER BY CLASS (top scoring, median):
//   Satellite             0.160   0.334   1.017   1.369
//   Unknown               0.099   0.276   0.821   1.017
//   Retroposon (SVA)      0.028   0.102   0.160   0.166
//   Simple_repeat         0.184   0.117   0.070   0.083
//   Low_complexity        0.061   0.047   0.054   0.074
//   LINE/SINE/LTR/DNA     0.000   0.007   0.028   0.051
//
// EFFECTIVENESS ANALYSIS:
//
// Most effective (score clearly separates repeat from non-repeat):
//   • STRs at 100–200bp — triTR_gt200 scores 2.53 at 100bp, diTR 1.20.
//     Signal *dilutes* at larger scales (triTR_gt200: 2.53→0.53 100→2000bp)
//     meaning multi-scale profiles correctly identify STR as local repeats.
//   • Tandem satellites at ≥1000bp — HSAT/CenSat/NewSat grow monotonically,
//     reaching 0.38–0.50 at 2000bp (vs 0.07 for non-repetitive).
//   • Telomeres at any scale — 1.90 at 100bp, 4.04 at 2000bp.
//   • RM/Unknown (tandem arrays) — strong monotonic growth to 1.02 at 2000bp.
//   • RM/Satellite — 0.16 at 100bp → 1.37 at 2000bp (true tandem arrays).
//
// Least effective (score provides little discrimination):
//   • Ancient interspersed repeats — LINE/SINE/LTR/DNA median 0.05–0.07 at
//     2000bp ≈ non-repetitive baseline. Alignment-based detection required.
//   • SegDups — 0.065 at 2000bp, indistinguishable from unique DNA.
//   • Short STRs in long windows — diTR_51to200: 1.20 at 100bp → 0.12 at 2000bp.
//
// INTERPRETATION GUIDE (per-k-mer score q):
//   q = 0      →  completely random / unique DNA
//   q < 0.1    →  no meaningful repeat structure (includes most RM elements)
//   0.1–0.6    →  weak to moderate repetition (STR-adjacent, diverged repeats)
//   0.6–1.0    →  moderately repetitive (short STRs, some satellites)
//   1.0–2.0    →  highly repetitive (long STRs, HSAT, compact satellites)
//   > 2.0      →  extremely repetitive (telomeres, long homopolymers)

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Calibration: multi-scale scores across all annotation sources",
          "[lancet][base][longdust_scorer][calibration]") {
  auto const ref_path = MakePath(FULL_DATA_DIR, CHM13_REF_NAME);
  REQUIRE(std::filesystem::exists(ref_path));
  auto const tsv_path = MakePath(FULL_DATA_DIR, CALIBRATION_TSV);
  REQUIRE(std::filesystem::exists(tsv_path));

  lancet::hts::Reference const ref(ref_path);
  lancet::base::LongdustQScorer const scorer(7, 10'001);
  constexpr double MARGIN = 0.001;

  auto const rows = LoadCalibrationTsv(tsv_path);
  INFO("Loaded " << rows.size() << " calibration rows from " << CALIBRATION_TSV);
  REQUIRE(!rows.empty());

  for (auto const& row : rows) {
    DYNAMIC_SECTION(row.mSource << " | " << row.mName << " @" << row.mScale << "bp") {
      auto const rgn = ref.MakeRegion(row.mRegion.c_str());
      CHECK(scorer.Score(rgn.SeqView()) == Catch::Approx(row.mExpectedScore).margin(MARGIN));
    }
  }
}
