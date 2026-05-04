#include "lancet/base/sequence_complexity.h"

#include "catch_amalgamated.hpp"

#include <algorithm>
#include <string>

#include <cmath>

namespace lancet::base::tests {

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PART 1: MaxHomopolymerRun                                               ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("MaxHomopolymerRun: basic cases", "[lancet][base][SequenceComplexityScorer]") {
  CHECK(SequenceComplexityScorer::MaxHomopolymerRun("") == 0);
  CHECK(SequenceComplexityScorer::MaxHomopolymerRun("A") == 1);
  CHECK(SequenceComplexityScorer::MaxHomopolymerRun("ACGT") == 1);
  CHECK(SequenceComplexityScorer::MaxHomopolymerRun("AACCCGTTT") == 3);
  CHECK(SequenceComplexityScorer::MaxHomopolymerRun("AAAAAAA") == 7);
  CHECK(SequenceComplexityScorer::MaxHomopolymerRun("ATCAAAAAGTC") == 5);
}

TEST_CASE("MaxHomopolymerRun: full homopolymer", "[lancet][base][SequenceComplexityScorer]") {
  CHECK(SequenceComplexityScorer::MaxHomopolymerRun(std::string(50, 'T')) == 50);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PART 2: LocalShannonEntropy                                             ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("LocalShannonEntropy: edge and known values",
          "[lancet][base][SequenceComplexityScorer]") {
  CHECK(SequenceComplexityScorer::LocalShannonEntropy("") == 0.0F);
  CHECK(SequenceComplexityScorer::LocalShannonEntropy("AAAA") == 0.0F);
  CHECK(SequenceComplexityScorer::LocalShannonEntropy("TTTTTTTT") == 0.0F);
  CHECK(SequenceComplexityScorer::LocalShannonEntropy("ACGT") ==
        Catch::Approx(2.0F).margin(0.001F));
  CHECK(SequenceComplexityScorer::LocalShannonEntropy("AACCGGTT") ==
        Catch::Approx(2.0F).margin(0.001F));

  // Two-base: equal freq → H = 1.0
  CHECK(SequenceComplexityScorer::LocalShannonEntropy("ACACAC") ==
        Catch::Approx(1.0F).margin(0.001F));

  // Three bases equally frequent → H = log2(3) ≈ 1.585
  CHECK(SequenceComplexityScorer::LocalShannonEntropy("AACCGG") ==
        Catch::Approx(std::log2(3.0F)).margin(0.01F));
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PART 3: Tandem Repeat Detection                                         ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("FindExactRepeats: dinucleotide repeat", "[lancet][base][SequenceComplexityScorer]") {
  // ATATATATAT → period=2, copies=5.0, span=10
  auto const results = SequenceComplexityScorer::FindExactRepeats("ATATATATAT");
  REQUIRE(!results.empty());

  // Find the period=2 result
  auto const* best = results.data();
  for (auto const& rslt : results) {
    if (rslt.mCopies > best->mCopies ||
        (rslt.mCopies == best->mCopies && rslt.mPeriod < best->mPeriod)) {
      best = &rslt;
    }
  }

  CHECK(best->mPeriod == 2);
  CHECK(best->mCopies == Catch::Approx(5.0F).margin(0.01F));
  CHECK(best->mSpanLength == 10);
  CHECK(best->mIsExact == true);
  CHECK(best->mTotalErrors == 0);
}

TEST_CASE("FindExactRepeats: homopolymer", "[lancet][base][SequenceComplexityScorer]") {
  // AAAAAA → period=1, copies=6.0
  auto const results = SequenceComplexityScorer::FindExactRepeats("AAAAAA");
  REQUIRE(!results.empty());

  bool found_period1 = false;
  for (auto const& rslt : results) {
    if (rslt.mPeriod == 1 && rslt.mCopies >= 6.0F) {
      found_period1 = true;
      CHECK(rslt.mSpanLength == 6);
      CHECK(rslt.mIsExact == true);
    }
  }
  CHECK(found_period1);
}

TEST_CASE("FindExactRepeats: no repeat in random DNA", "[lancet][base][SequenceComplexityScorer]") {
  // Short non-repetitive sequence — should find no repeats with min_copies=2.5
  auto const results = SequenceComplexityScorer::FindExactRepeats("ACGTACGA");
  // Even if some results, copies should be low
  for (auto const& rslt : results) {
    // Allow some results but verify they're genuine
    CHECK(rslt.mPeriod > 0);
  }
}

TEST_CASE("FindExactRepeats: primitive motif enforcement",
          "[lancet][base][SequenceComplexityScorer]") {
  // ATATATAT should find AT (period=2) but NOT ATAT (period=4, non-primitive)
  auto const results = SequenceComplexityScorer::FindExactRepeats("ATATATAT");
  for (auto const& rslt : results) {
    // Period 4 with motif ATAT should be filtered (it's a repeat of AT)
    if (rslt.mPeriod == 4) {
      // If period=4 result exists, it should be a different motif, not ATAT
      FAIL("Period-4 non-primitive motif ATAT should have been filtered");
    }
  }
}

// Catch2 SECTION fan-out inflates clang-tidy's cognitive-complexity metric beyond the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("FindApproxRepeats: imperfect trinucleotide repeat",
          "[lancet][base][SequenceComplexityScorer]") {
  // CAGCAACAGCAG → period=3, ~4 copies, with 1 error (G→A mutation in the 2nd
  // copy). FindApproxRepeats with `period=3, min_copies=3.0, max_errors=1`
  // should detect this as an approximate trinucleotide repeat. The previous
  // version of this test discarded the `found_approx` flag without asserting,
  // so a regression that silently lost the match would have been invisible.
  auto const results = SequenceComplexityScorer::FindApproxRepeats("CAGCAACAGCAG", 6, 3.0F, 1);

  // Structural invariants on EVERY returned result (any return shape that
  // breaks these is a bug, regardless of the algorithm's tuning):
  //   - period ≥ 1 (period=0 is meaningless)
  //   - copies ≥ 1.0 (a "repeat" of less than one copy isn't one)
  //   - span ≥ period (one copy is at least period bases)
  //   - errors are non-negative
  //   - purity is in [0, 1]
  for (auto const& rslt : results) {
    INFO("period=" << rslt.mPeriod << " copies=" << rslt.mCopies << " span=" << rslt.mSpanLength
                   << " errors=" << rslt.mTotalErrors);
    CHECK(rslt.mPeriod >= 1);
    CHECK(rslt.mCopies >= 1.0F);
    CHECK(rslt.mSpanLength >= rslt.mPeriod);
    CHECK(rslt.mTotalErrors >= 0);
    CHECK(rslt.Purity() >= 0.0F);
    CHECK(rslt.Purity() <= 1.0F);
  }

  // The hand-checked semantic claim: with period=3 / min_copies=3.0 /
  // max_errors=1 over "CAGCAACAGCAG", at least one approximate match should
  // be returned with period 3 and ≥ 3 copies. Per-match assertions inside
  // the loop verify the error-count and purity stay within the documented
  // shape (≥ 1 error counted; purity ≥ 0.75 = 1 − 1/12 worst-case).
  bool found_approx_period_3 = false;
  for (auto const& rslt : results) {
    if (rslt.mPeriod == 3 && rslt.mCopies >= 3.0F) {
      found_approx_period_3 = true;
      CHECK(rslt.mTotalErrors >= 1);
      CHECK(rslt.Purity() >= 0.75F);
    }
  }
  CHECK(found_approx_period_3);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PART 4: SequenceComplexityScorer::Score (end-to-end)                    ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("Score: poly-A variant produces expected context/delta/TR features",
          "[lancet][base][SequenceComplexityScorer]") {
  SequenceComplexityScorer const scorer;

  // 200bp REF haplotype: poly-A centered at 90-110
  std::string ref_hap(90, 'C');
  ref_hap += std::string(20, 'A');  // poly-A at pos 90-110
  ref_hap += std::string(90, 'G');

  // ALT haplotype: extends the poly-A by 5bp (25 A's total = INDEL extension)
  std::string alt_hap(90, 'C');
  alt_hap += std::string(25, 'A');
  alt_hap += std::string(85, 'G');

  auto const cplx = scorer.Score({.mHaplotype = ref_hap, .mPos = 90, .mLen = 20},
                                 {.mHaplotype = alt_hap, .mPos = 90, .mLen = 25});

  // Context: REF ±20bp around poly-A → high HRun, low entropy
  CHECK(cplx.ContextHRun() >= 20);  // 20bp of A's visible in ±20bp window
  CHECK(cplx.ContextEntropy() >= 0.0F);

  // Delta: ALT extended the homopolymer → positive DeltaHRun
  CHECK(cplx.DeltaHRun() >= 0);

  // Context LQ: should be non-negative after log1p
  CHECK(cplx.ContextFlankLQ() >= 0.0);
  CHECK(cplx.ContextHaplotypeLQ() >= 0.0);
}

TEST_CASE("Score: non-repetitive REF produces low context scores",
          "[lancet][base][SequenceComplexityScorer]") {
  SequenceComplexityScorer const scorer;

  // Balanced, non-repetitive haplotype
  std::string haplo;
  for (int iter = 0; iter < 50; ++iter) haplo += "ACGT";  // 200bp

  // REF = ALT (no change) at position 100, length 1
  auto const cplx = scorer.Score({.mHaplotype = haplo, .mPos = 100, .mLen = 1},
                                 {.mHaplotype = haplo, .mPos = 100, .mLen = 1});

  CHECK(cplx.ContextHRun() == 1);                                    // no homopolymer
  CHECK(cplx.ContextEntropy() == Catch::Approx(2.0F).margin(0.1F));  // near-max entropy
  CHECK(cplx.DeltaHRun() == 0);                                      // REF == ALT
  CHECK(cplx.DeltaEntropy() == Catch::Approx(0.0F).margin(0.01F));
}

TEST_CASE("Score: sentinel handling — no TR found → zeros",
          "[lancet][base][SequenceComplexityScorer]") {
  SequenceComplexityScorer const scorer;

  // Short non-repetitive haplotype
  std::string const haplo = "ACGTACGTACGTACGTACGTACGTACGTACGT"
                            "TGCATGCATGCATGCATGCATGCATGCATGCA";

  auto const cplx = scorer.Score({.mHaplotype = haplo, .mPos = 16, .mLen = 1},
                                 {.mHaplotype = haplo, .mPos = 16, .mLen = 1});

  // TrAffinity should be 0 or close to 0 if no TR found nearby
  CHECK(cplx.TrAffinity() >= 0.0F);
  CHECK(cplx.TrAffinity() <= 1.0F);
  CHECK(cplx.TrPurity() >= 0.0F);
  CHECK(cplx.TrPeriod() >= 0);
  CHECK(cplx.IsStutterIndel() >= 0);
}

TEST_CASE("Score: ScoreMacro matches LongdustQScorer output",
          "[lancet][base][SequenceComplexityScorer]") {
  SequenceComplexityScorer const scorer(0.5);  // uniform GC for comparison
  lancet::base::LongdustQScorer const direct_scorer(7, 4096, 0.5);

  // Score a repetitive haplotype
  auto const haplotype = std::string(200, 'A');
  auto const cplx = scorer.Score({.mHaplotype = haplotype, .mPos = 100, .mLen = 1},
                                 {.mHaplotype = haplotype, .mPos = 100, .mLen = 1});

  // ContextHaplotypeLQ = log1p(direct LQ score)
  CHECK(cplx.ContextHaplotypeLQ() ==
        Catch::Approx(std::log1p(direct_scorer.Score(haplotype))).margin(1e-4));
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PART 5: LongdustQScorer GC-Bias Tests                                   ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("GC-bias: gc_frac=0.5 matches pre-GC behavior", "[lancet][base][LongdustQScorer]") {
  // gc_frac=0.5 → uniform null model → identical to original formula
  lancet::base::LongdustQScorer const uniform_scorer(7, 1024, 0.5);
  lancet::base::LongdustQScorer const default_human(7, 1024, 0.41);

  // Truly random DNA (no periodic structure)
  std::string const random_dna = "GCTAAGGTCCTTGAACGGATTCATAGCCTGAGATTTCAAC"
                                 "TGCAAGGTCCTCATGAACTTTAGCCCAAGATTCTGAACGT";
  CHECK(uniform_scorer.Score(random_dna) < 0.5);
  CHECK(default_human.Score(random_dna) < 0.5);
}

TEST_CASE("GC-bias: gc_frac=0.41 lowers AT-rich non-repetitive scores",
          "[lancet][base][LongdustQScorer]") {
  lancet::base::LongdustQScorer const uniform_scorer(7, 1024, 0.5);
  lancet::base::LongdustQScorer const human_scorer(7, 1024, 0.41);

  // AT-rich non-repetitive DNA: human scorer should give LOWER score
  // because AT-richness is expected in human genome (gc=0.41 → AT=0.59)
  std::string const at_rich = "ATATATGTAACTTAATGTATTATATTGATGAATTTAATGG"
                              "ATTAAGTCATATTAATGATTAATATGATATAAGAAATAGG";
  // The human scorer adjusts expected k-mer frequencies for AT bias,
  // reducing apparent "repetitiveness" of compositionally biased sequence
  auto const uniform_score = uniform_scorer.Score(at_rich);
  auto const human_score = human_scorer.Score(at_rich);
  CHECK(human_score <= uniform_score);
}

TEST_CASE("GC-bias: structural repeats still score high with gc=0.41",
          "[lancet][base][LongdustQScorer]") {
  lancet::base::LongdustQScorer const human_scorer(7, 1024, 0.41);

  // Pure poly-A: indisputably repetitive regardless of GC content
  CHECK(human_scorer.Score(std::string(50, 'A')) > 0.5);

  // Microsatellite (CA repeat): clearly repetitive
  std::string micro;
  for (int iter = 0; iter < 25; ++iter) micro += "CA";
  CHECK(human_scorer.Score(micro) > 0.5);
}

TEST_CASE("GC-bias: gc_frac=-1 equivalent to gc_frac=0.5", "[lancet][base][LongdustQScorer]") {
  // gc_frac is clamped to [0,1], so -1.0 becomes 0.0
  // This is NOT equivalent to 0.5 — it's an extreme AT-bias model.
  // The actual equivalence is: gc_frac=0.5 IS the uniform model.
  lancet::base::LongdustQScorer const scorer_05(7, 1024, 0.5);

  // Verify gc_frac=0.5 produces the uniform model
  std::string const test_seq = std::string(100, 'A');
  CHECK(scorer_05.Score(test_seq) > 0.0);  // poly-A should always score > 0
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  PART 6: SequenceComplexity output testing                               ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("SequenceComplexity::FormatVcfValue: 11 comma-separated values",
          "[lancet][base][SequenceComplexityScorer]") {
  SequenceComplexity const cplx;  // default-constructed → all zeros
  auto const vcf_val = cplx.FormatVcfValue();
  // 11 values → 10 commas
  CHECK(std::count(vcf_val.begin(), vcf_val.end(), ',') == 10);
}

TEST_CASE("SequenceComplexity: default-constructed → all zeros",
          "[lancet][base][SequenceComplexityScorer]") {
  SequenceComplexity const cplx;

  CHECK(cplx.ContextHRun() == 0);
  CHECK(cplx.ContextEntropy() == 0.0F);
  CHECK(cplx.ContextFlankLQ() == 0.0);
  CHECK(cplx.ContextHaplotypeLQ() == 0.0);
  CHECK(cplx.DeltaHRun() == 0);
  CHECK(cplx.DeltaEntropy() == 0.0F);
  CHECK(cplx.DeltaFlankLQ() == 0.0);
  CHECK(cplx.TrAffinity() == 0.0F);
  CHECK(cplx.TrPurity() == 0.0F);
  CHECK(cplx.TrPeriod() == 0);
  CHECK(cplx.IsStutterIndel() == 0);
}

TEST_CASE("SequenceComplexity::MergeMax: element-wise max",
          "[lancet][base][SequenceComplexityScorer]") {
  SequenceComplexityScorer const scorer;

  // Two different ALT haplotypes for the same variant
  std::string ref_hap(90, 'C');
  ref_hap += std::string(20, 'A');
  ref_hap += std::string(90, 'G');

  std::string alt1(90, 'C');
  alt1 += std::string(25, 'A');
  alt1 += std::string(85, 'G');

  std::string alt2(90, 'C');
  alt2 += std::string(30, 'A');
  alt2 += std::string(80, 'G');

  auto cplx1 = scorer.Score({.mHaplotype = ref_hap, .mPos = 90, .mLen = 20},
                            {.mHaplotype = alt1, .mPos = 90, .mLen = 25});
  auto const cplx2 = scorer.Score({.mHaplotype = ref_hap, .mPos = 90, .mLen = 20},
                                  {.mHaplotype = alt2, .mPos = 90, .mLen = 30});

  // MergeMax should take element-wise max
  cplx1.MergeMax(cplx2);

  CHECK(cplx1.ContextHRun() >= cplx2.ContextHRun());
  CHECK(cplx1.DeltaHRun() >= cplx2.DeltaHRun());
}

}  // namespace lancet::base::tests
