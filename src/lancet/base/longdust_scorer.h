#ifndef SRC_LANCET_BASE_LONGDUST_SCORER_H_
#define SRC_LANCET_BASE_LONGDUST_SCORER_H_

#include "lancet/base/rev_comp.h"
#include "lancet/base/types.h"

#include "absl/strings/str_cat.h"
#include "spdlog/fmt/bundled/core.h"

#include <algorithm>
#include <array>
#include <numbers>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cmath>
#include <cstdint>

namespace lancet::base {

// ============================================================================
// LongdustQScorer — k-mer complexity scorer for DNA sequences
//
// PROVENANCE
// ----------
// This scorer implements the Q(x) complexity formula from:
//
//   Li, H. "Finding low-complexity filter with longdust" (2025)
//   https://arxiv.org/abs/2509.07357  |  https://github.com/lh3/longdust
//
// Longdust is a tool for finding low-complexity *intervals* in whole genomes.
// It uses the Q(x) formula internally as part of a sliding-window scanner
// with X-drop extension, minimum-start-count triggers, and approximation
// heuristics. We extract just the scoring formula and apply it directly to
// known subsequences around variants. The formula, f(ℓ) precomputation,
// k-mer encoding, and both-strand logic are verified to produce identical
// results via cross-validation against longdust's C implementation
// (see tests/base/longdust_scorer_test.cpp).
//
// GC-BIAS CORRECTION
// ------------------
// The scorer supports an optional GC-bias correction via the gc_frac
// constructor parameter. By default, gc_frac=0.41 (human genome-wide
// average GC content, confirmed by Lander et al. 2001, Piovesan et al.
// 2019, and Nurk et al. 2022 T2T-CHM13).
//
// Without correction (gc_frac=0.5), the null model assumes uniform base
// frequencies (25% each). In AT-rich or GC-rich regions, this inflates
// Q(x) for non-repetitive but compositionally biased DNA. The correction
// groups 4^k possible k-mers into k+1 binomial equivalence classes based
// on their GC content, adjusting the expected count λ for each class.
// This only changes f(ℓ) precomputation — the hot-path Score() method
// remains completely unchanged (zero runtime overhead).
//
// IMPORTANT: gc_frac must be a GLOBAL or broad regional background
// fraction, NOT the local GC of the scored window. A poly-A insertion
// would locally compute as 0% GC, causing the scorer to expect all-A
// k-mers and produce Q(x)≈0, making it blind to the repeat.
//
// Set gc_frac=0.5 to disable GC correction (uniform model).
// For non-human genomes, set to the organism's genome-wide GC fraction
// (e.g., ~0.20 for P. falciparum, ~0.65 for S. cerevisiae).
//
// ============================================================================
//
// INTUITION: WHAT IS COMPLEXITY SCORING?
// ======================================
//
// DNA complexity scoring answers: "how repetitive is this sequence?"
//
// The key idea is to count k-mers (short words of length k) and see if any
// k-mer appears unusually often. In random DNA, each k-mer appears roughly
// the same number of times. In repetitive DNA, a few k-mers dominate:
//
//   Random DNA (100bp, k=7):
//     ATCGATCGTAGCTGATCGATCGTAGC...
//     → 94 k-mers, each appearing 0 or 1 times
//     → score ≈ 0  (uniform distribution, no concentration)
//
//   Homopolymer run (100bp, k=7):
//     AAAAAAAAAAAAAAAAAAAAAAAAA...
//     → only one k-mer (AAAAAAA) appears 94 times
//     → score >> 0  (extreme concentration into one k-mer)
//
//   Dinucleotide repeat (100bp, k=7):
//     CACACACACACACACACACACACACA...
//     → only two k-mers (CACACAC, ACACACA) each appear ~47 times
//     → score > 0  (high concentration into two k-mers)
//
// The score quantifies this "concentration" by comparing the observed
// k-mer count distribution against a random (Poisson) null model.
//
// ============================================================================
//
// THE MATH: Q(x) COMPLEXITY SCORE
// ================================
//
// Given a DNA string x, let:
//   k       = k-mer size (we use k=7, giving 4^7 = 16,384 possible k-mers)
//   c_x(t)  = how many times k-mer t appears in x
//   ℓ       = total number of valid k-mers extracted from x
//   4^k     = total possible distinct k-mers
//
// Step 1: Measure k-mer concentration
// ------------------------------------
// Compute the sum of log-factorials of all k-mer counts:
//
//   Σ_t log(c_x(t)!)
//
// Why log-factorials? If N items are distributed among B bins:
//   - Perfectly uniform (each bin has N/B):  Σ log((N/B)!) is small
//   - All in one bin (one has N, rest 0):    log(N!) is large
//
// Example (ℓ=10 k-mers, k=7):
//   Random:     counts = [1,1,1,1,1,1,1,1,1,1, 0,0,...,0]
//               Σ log(c!) = 10 × log(1!) = 0
//
//   Homopolymer: counts = [10, 0,0,...,0]
//               Σ log(c!) = log(10!) = 15.1
//
// Step 2: Subtract the expected value under the null model
// ---------------------------------------------------------
// Under a random (Poisson) model, each k-mer has expected count λ = ℓ/4^k.
// The expected value of Σ log(c!) under this model is f(ℓ):
//
//   f(ℓ) = 4^k × f_single(ℓ / 4^k)
//
//   f_single(λ) = E[log(N!)] where N ~ Poisson(λ)
//               = e^{-λ} × Σ_{n=2}^∞ log(n!) × λ^n / n!
//
// For large λ (≥ 30), Stirling's approximation is used.
//
// The raw complexity score is:
//
//   ╔═══════════════════════════════════════╗
//   ║  Q(x) = Σ_t log(c_x(t)!) − f(ℓ)    ║
//   ╚═══════════════════════════════════════╝
//
// Q(x) = 0 means the k-mer distribution matches the random expectation.
// Q(x) > 0 means k-mers are more concentrated than expected (repetitive).
//
// Step 3: Normalize
// ------------------
// We report the per-k-mer score:
//
//   q(x) = max(0, Q(x) / ℓ)
//
// This makes scores comparable across different sequence lengths.
//
// Interpretation:
//   q ≈ 0    →  random / unique sequence
//   q > 0.6  →  moderately repetitive (longdust's default LCR threshold)
//   q > 1.0  →  highly repetitive (strong tandem repeat signal)
//   q > 2.0  →  extremely repetitive (long homopolymer, satellite DNA)
//
// ============================================================================
//
// BOTH-STRAND SCORING
// ===================
//
// Repeats can be strand-asymmetric in k-mer space. For example, a poly-T
// run uses k-mer TTTTTTT on the forward strand, but its reverse complement
// (poly-A) uses AAAAAAA — a completely different k-mer. To handle this,
// we score both the forward and reverse complement, returning the maximum:
//
//   Score(x) = max( q(x_forward), q(x_revcomp) )
//
// This matches longdust's ld_dust2 behavior.
//
// ============================================================================
//
// MULTI-SCALE VARIANT ANNOTATION
// ==============================
//
// Each variant is scored at 5 flanking distances to capture repeats at
// different biological scales:
//
//   Haplotype: [==========|VVVV|==========]
//                         ^    ^
//                      var_pos  var_pos + var_len
//
//   Scale 0 (5bp):    [---5--|VVVV|---5--]   → homopolymer/short repeat
//   Scale 1 (10bp):   [--10--|VVVV|--10--]   → homopolymer/dinucleotide
//   Scale 2 (50bp):   [--50--|VVVV|--50--]   → microsatellite/STR
//   Scale 3 (100bp):  [-100--|VVVV|-100--]   → longer STR/VNTR
//   Scale 4 (full):   [entire haplotype   ]  → satellite/assembly quality
//
//   Clamped: start = max(0, var_pos − flank)
//            end   = min(hap_len, var_pos + var_len + flank)
//
//   VCF output: ALT_LCR=q5,q10,q50,q100,qfull
//   Example:    ALT_LCR=2.1,1.8,0.7,0.3,0.1
//
// ============================================================================
//
// EMPIRICAL SCORE BEHAVIOUR (from CHM13v2.0 calibration, k=7)
// ============================================================
//
// Balanced analysis over 133,714 regions (≤5000 per subcategory, evenly spaced)
// from GIAB stratifications, RepeatMasker (17 classes), CenSat, Telomere,
// NewSatellite, and CompositeRepeats annotations.
// Scales: flank = [50, 100, 500, 1000] → window = [100, 200, 1000, 2000] bp
//
// Full data in:
//   data/lcr_calibration_balanced_scored.tsv.gz  (534,856 rows)
//   data/lcr_analysis_by_category.txt            (category breakdown)
//   data/lcr_analysis_by_subcategory.txt         (subcategory breakdown)
//
// CATEGORY-LEVEL SUMMARY (median scores)
// ──────────────────────────────────────
//
//   Source Category       100bp   200bp  1000bp  2000bp
//   ─────────────────── ─────── ─────── ─────── ───────
//   Telomere (TTAGGG)n    1.897   2.554   3.678   4.037
//   NewSatellite          0.138   0.286   0.273   0.498
//   GIAB/LowComplexity    0.236   0.138   0.081   0.094
//   CenSat                0.005   0.017   0.113   0.174
//   RepeatMasker (all)    0.005   0.017   0.052   0.079
//   GIAB/Union            0.000   0.010   0.048   0.073
//   CompositeRepeats      0.005   0.010   0.051   0.069
//   SegDups               0.005   0.010   0.039   0.065
//   GIAB/OtherDifficult   0.005   0.007   0.182   0.331
//
// KEY GIAB SUBCATEGORIES (median scores, sorted by 2000bp desc)
// ─────────────────────────────────────────────────────────────
//   triTR_gt200             2.527   3.212   0.929   0.529
//   quadTR_gt200            2.262   2.937   0.755   0.386
//   satellites              0.013   0.037   0.240   0.375
//   triTR_51to200           1.295   0.785   0.288   0.198
//   quadTR_51to200          1.198   0.804   0.245   0.179
//   diTR_51to200            1.198   0.595   0.155   0.121
//   quadTR_20to50           0.259   0.142   0.085   0.095
//   homopolymer_gt20        0.399   0.201   0.088   0.094
//   diTR_11to50             0.172   0.131   0.076   0.089
//   homopolymer_gt11        0.184   0.097   0.068   0.083
//   triTR_15to50            0.126   0.073   0.060   0.080
//   homopolymer_7to11       0.013   0.017   0.047   0.072
//   homopolymer_4to6        0.000   0.010   0.038   0.064
//   SegDups                 0.005   0.010   0.039   0.065
//
// REPEATMASKER BY CLASS (median scores, sorted by 2000bp desc)
// ────────────────────────────────────────────────────────────
//   Satellite             0.160   0.334   1.017   1.369
//   Unknown               0.099   0.276   0.821   1.017
//   Retroposon (SVA)      0.028   0.102   0.160   0.166
//   Beta                  0.000   0.014   0.057   0.097
//   Simple_repeat         0.184   0.117   0.070   0.083
//   Low_complexity        0.061   0.047   0.054   0.074
//   SINE                  0.000   0.007   0.049   0.072
//   LINE                  0.000   0.010   0.041   0.067
//   DNA                   0.000   0.007   0.033   0.057
//   LTR                   0.000   0.007   0.028   0.051
//
// WHERE THE SCORE IS MOST EFFECTIVE
// ─────────────────────────────────
// The score discriminates best when the window captures ≥3 copies of the
// repeat unit, allowing concentration to exceed the Poisson null:
//
//   • Di/tri/quad-nucleotide STRs → strong at 100–200bp, then dilute as
//     window grows (triTR_gt200: 2.53 at 100bp → 0.53 at 2000bp). The
//     score correctly captures that STR signal is local.
//   • Tandem satellites (HSAT/CenSat) → emerge at 1000bp, peak at 2000bp
//   • Telomeres → unmistakable at all scales (1.90 at 100bp, 4.04 at 2000bp)
//   • RM/Unknown, RM/Satellite → grow monotonically with scale (1.02, 1.37
//     at 2000bp). These are tandem satellite arrays correctly detected.
//   • Homopolymers → detectable at 100bp (homopolymer_gt20 = 0.40)
//
// WHERE THE SCORE IS LEAST EFFECTIVE
// ──────────────────────────────────
//   • Ancient interspersed repeats (LINE, SINE, LTR): single-copy, diverged
//     insertions. At 2000bp, median scores are 0.051–0.072 — barely above
//     the non-repetitive baseline. Detection requires alignment, not k-mer
//     concentration.
//
//   • Segmental duplications: median 0.065 at 2000bp, indistinguishable from
//     non-repetitive DNA. SegDup detection needs whole-genome alignment.
//
//   • Short STRs in long windows: a diTR_51to200 scores 1.20 at 100bp but
//     only 0.12 at 2000bp — the signal is diluted by flanking unique
//     sequence. Multi-scale scoring captures this correctly.
//
// SCORE DISTRIBUTION (balanced dataset, 133,714 regions per scale)
// ────────────────────────────────────────────────────────────────
//    100bp: 73% nonzero, 40% >0.1, 12% >0.6,  8% >1.0
//    200bp: 88% nonzero, 37% >0.1,  8% >0.6,  3% >1.0
//   1000bp: ~100% nonzero, 33% >0.1,  8% >0.6,  6% >1.0
//   2000bp: 100% nonzero, 41% >0.1,  8% >0.6,  6% >1.0
//
// SCALE SELECTION GUIDANCE
// ───────────────────────
//   • STR/microsatellite filtering: 100bp scale is most discriminating
//   • Satellite/centromeric filtering: 1000–2000bp scales needed
//   • The multi-scale profile is a fingerprint: STRs peak at small scales
//     and decay, while satellites grow monotonically with window size
//
// ============================================================================
//
// COMPARISON WITH LONGDUST
// ========================
//
// What is identical:
//   - k-mer size (k=7) and base encoding (A=0, C=1, G=2, T=3)
//   - The Q(x) formula: Σ log(c!) − f(ℓ)
//   - f(ℓ) computation: Poisson series + Stirling approximation for λ≥30
//   - Both-strand scoring: max(forward, reverse complement)
//   - GC correction support: longdust uses gc=-1 (uniform) by default.
//     LongdustQScorer supports an optional gc_frac parameter for GC-bias
//     correction via a binomial null model (see GC-BIAS CORRECTION above).
//
// What differs:
//   - Longdust scans long sequences with a sliding window (ws=5000) to
//     find LCR *intervals* using backward/forward sweeps with X-drop
//     extension (xdrop_len=50), minimum start count (min_start_cnt=3),
//     and optional O(Lw) approximation. It outputs intervals, not scores.
//   - LongdustQScorer computes Q(x)/ℓ directly on a known subsequence (the
//     flanking region around a variant). No scanning is needed because
//     the variant location is already known from graph assembly + MSA.
//   - LongdustQScorer outputs the raw per-k-mer score, not intervals.
//
// ============================================================================

// NOTE: Scale configuration constants (NUM_LCR_SCALES, LCR_FLANKS, etc.) have
// been moved to the SequenceComplexity class which manages multi-scale scoring.
// The LongdustQScorer is now a pure scoring primitive — it takes a sequence and
// returns the Q(x) complexity score. Scale management is handled by the caller.

/// k-mer sizes used by SequenceComplexity for different scales.
static constexpr int LONGDUST_FLANK_K = 4;
static constexpr int LONGDUST_HAPLOTYPE_K = 7;

// ============================================================================
// FormatComplexityScore: format a floating-point score for VCF output.
//
//   - Up to 3 decimal places, no trailing zeros
//   - Examples: 0.0→"0", 1.5→"1.5", 0.123→"0.123", 2.100→"2.1"
// ============================================================================
inline auto FormatComplexityScore(f64 val) -> std::string {
  auto txt = fmt::format("{:.3f}", val);
  if (txt.find('.') != std::string::npos) {
    auto const last_nonzero = txt.find_last_not_of('0');
    txt.erase(last_nonzero + 1);
    if (txt.back() == '.') {
      txt.pop_back();
    }
  }
  return txt;
}

/// Format a span of scores as a comma-separated VCF INFO value.
template <typename Container>
inline auto FormatComplexityScores(Container const& scores) -> std::string {
  std::string result;
  bool first = true;
  for (auto const& score : scores) {
    if (!first) {
      absl::StrAppend(&result, ",");
    }
    absl::StrAppend(&result, FormatComplexityScore(static_cast<f64>(score)));
    first = false;
  }
  return result;
}

// ============================================================================
// Base encoding table
//
// Maps each ASCII character to a 2-bit DNA code:
//   A/a → 0 (binary 00)
//   C/c → 1 (binary 01)
//   G/g → 2 (binary 10)
//   T/t → 3 (binary 11)
//   Everything else → 4 (sentinel: breaks the k-mer window)
//
// This is packed into a k-mer integer by shifting and OR-ing:
//   Example (k=3): sequence ACG → 00·01·10 = binary 000110 = decimal 6
// ============================================================================
constexpr auto MakeDnaEncodeTable() -> std::array<u8, 256> {
  std::array<u8, 256> tbl{};
  for (auto& val : tbl) {
    val = 4;
  }
  tbl['A'] = 0;
  tbl['a'] = 0;
  tbl['C'] = 1;
  tbl['c'] = 1;
  tbl['G'] = 2;
  tbl['g'] = 2;
  tbl['T'] = 3;
  tbl['t'] = 3;
  return tbl;
}

inline constexpr std::array<u8, 256> DNA_ENCODE_TABLE = MakeDnaEncodeTable();

class LongdustQScorer {
 public:
  /// @param k       k-mer size (default 7: 4^7 = 16,384 possible k-mers)
  /// @param max_len Maximum sequence length for precomputing f(ℓ)
  /// @param gc_frac Global background GC fraction for compositional bias correction.
  ///                Default 0.41 (human genome-wide average: Lander et al. 2001,
  ///                Piovesan et al. 2019, Nurk et al. 2022). Set to 0.5 for
  ///                uniform (no correction). For non-human genomes, set to
  ///                the organism's genome-wide GC (e.g., 0.20 for P. falciparum).
  explicit LongdustQScorer(int k = 7, int max_len = 1024, f64 gc_frac = 0.41)
      : mK(k),
        mMask((1U << (2 * k)) - 1),
        mNumKmers(1U << (2 * k)),
        mGc(std::clamp(gc_frac, 0.0, 1.0)) {
    PrecomputeF(max_len);
  }

  // ============================================================================
  // Score: both-strand complexity score
  //
  // Computes q(x) on both the forward sequence and its reverse complement,
  // returns the maximum. This ensures strand-symmetric detection of repeats.
  //
  // Example:
  //   Forward:  TTTTTTTTTTTTT  →  k-mer {TTTTTTT} count=7  →  q=1.2
  //   RevComp:  AAAAAAAAAAAAA  →  k-mer {AAAAAAA} count=7  →  q=1.2
  //   Score = max(1.2, 1.2) = 1.2
  //
  //   But for an asymmetric repeat like ATATATATATAT:
  //   Forward k-mers:  {ATATATA, TATATAT} each ~3-4x  →  q_fwd
  //   RevComp k-mers:  {ATATATA, TATATAT} each ~3-4x  →  q_rev
  //   Score = max(q_fwd, q_rev)
  // ============================================================================
  [[nodiscard]] auto Score(std::string_view seq) const -> f64 {
    f64 const fwd_score = ScoreOneStrand(seq);
    f64 const rev_score = ScoreOneStrand(RevComp(seq));
    return std::max(fwd_score, rev_score);
  }

  // ============================================================================
  // ScoreOneStrand: single-strand complexity score q(x) = Q(x) / ℓ
  //
  // Walks through the sequence in four steps:
  //
  //   Step 1: Count k-mers
  //   ─────────────────────
  //   Slide a window of width k across the sequence. Each base is encoded
  //   as 2 bits and shifted into a k-mer integer. N bases reset the window.
  //
  //     seq:    A  T  C  G  A  T  C  G  A  T  C  (length=11, k=7)
  //             └────k-mer 1────┘
  //                └────k-mer 2────┘
  //                   └────k-mer 3────┘
  //                      └────k-mer 4────┘
  //                         └────k-mer 5────┘
  //                                          ℓ = 5 valid k-mers
  //
  //   Step 2: Sum log-factorials
  //   ──────────────────────────
  //   For each k-mer with count ≥ 2, add log(count!).
  //
  //     Example (homopolymer, ℓ=10):
  //       Only k-mer AAAAAAA has count=10.
  //       Σ log(c!) = log(10!) = 15.1
  //
  //     Example (random DNA, ℓ=94):
  //       Each k-mer appears 0 or 1 times.
  //       Σ log(c!) = 0  (log(0!)=log(1!)=0)
  //
  //   Step 3: Subtract f(ℓ)
  //   ──────────────────────
  //   f(ℓ) is the expected value of Σ log(c!) under the Poisson null model.
  //   Q(x) = Σ log(c!) − f(ℓ)
  //
  //   Step 4: Normalize
  //   ─────────────────
  //   q(x) = max(0, Q(x) / ℓ)
  //
  // Complexity: O(|seq|) time, O(4^k) space (32KB for k=7, fits L1 cache).
  // ============================================================================
  [[nodiscard]] auto ScoreOneStrand(std::string_view seq) const -> f64 {
    auto const n_positions = static_cast<i64>(seq.size()) - mK + 1;
    if (n_positions <= 0) {
      return 0.0;
    }

    // Step 1: Count k-mers in a flat array indexed by the 2k-bit k-mer code
    thread_local std::vector<u16> counts;
    counts.assign(mNumKmers, 0);

    u32 kmer = 0;
    int run = 0;
    usize valid_kmers = 0;

    for (char const chr : seq) {
      auto const base = DNA_ENCODE_TABLE[static_cast<u8>(chr)];
      if (base < 4) {
        kmer = ((kmer << 2) | base) & mMask;
        if (++run >= mK) {
          counts[kmer]++;
          valid_kmers++;
        }
      } else {
        run = 0;  // N or invalid base resets the k-mer window
      }
    }

    if (valid_kmers == 0) {
      return 0.0;
    }

    // Step 2: Σ_t log(c_x(t)!) — only counts ≥ 2 contribute (log(0!)=log(1!)=0)
    f64 sum_log_factorial = 0.0;
    for (u32 idx = 0; idx < mNumKmers; ++idx) {
      if (counts[idx] >= 2) {
        sum_log_factorial += std::lgamma(static_cast<f64>(counts[idx] + 1));
      }
    }

    // Step 3: Q(x) = Σ log(c!) − f(ℓ)
    auto const ell = static_cast<int>(valid_kmers);
    f64 const f_val = (std::cmp_less(ell, mF.size())) ? mF[ell] : ComputeF(ell);
    f64 const q_score = sum_log_factorial - f_val;

    // Step 4: normalize, clamp at 0
    return std::max(0.0, q_score / static_cast<f64>(valid_kmers));
  }

 private:
  int mK;               // k-mer size (e.g., 7 → 7-mers)
  u32 mMask;            // bitmask to extract a k-mer from the rolling integer: (1 << 2k) - 1
  u32 mNumKmers;        // total possible k-mers: 4^k (e.g., 16,384 for k=7)
  f64 mGc;              // global background GC fraction for bias correction
  std::vector<f64> mF;  // precomputed f(ℓ) table: mF[ℓ] = expected Σ log(c!) under null model

  // ============================================================================
  // ComputeFSingle: expected log-factorial per k-mer under Poisson(λ)
  //
  // Computes f_single(λ) = E[log(N!)] where N ~ Poisson(λ).
  //
  // This is the expected "concentration signal" for a single k-mer when
  // ℓ k-mers are distributed randomly among 4^k bins (each bin gets
  // λ = ℓ/4^k expected count).
  //
  // Two regimes:
  //   λ < 30:  Exact series — sum the Poisson-weighted log-factorials:
  //            f_single(λ) = e^{-λ} × Σ_{n=2}^∞ log(n!) × λ^n/n!
  //            Converges quickly since Poisson PMF decays super-exponentially.
  //
  //   λ ≥ 30:  Stirling's approximation (faster, still accurate):
  //            f_single(λ) ≈ 0.5·log(2πeλ) − 1/(12λ)·(1 + 0.5/λ + 19/(30λ²))
  //                          + λ·(log(λ) − 1)
  //
  // Example values:
  //   λ = 0.01:  f_single ≈ 0.0   (most k-mers appear 0 or 1 times)
  //   λ = 1.0:   f_single ≈ 0.19  (some k-mers appear 2+, modest signal)
  //   λ = 10.0:  f_single ≈ 14.0  (significant expected concentration)
  // ============================================================================
  [[nodiscard]] static auto ComputeFSingle(f64 lambda) -> f64 {
    if (lambda < 1e-10) {
      return 0.0;
    }

    if (lambda >= 30.0) {
      // Stirling's approximation
      f64 const inv_lambda = 1.0 / lambda;
      f64 const stirling_val =
          (0.5 * std::log(2.0 * std::numbers::pi * std::numbers::e * lambda)) -
          (inv_lambda /
           12.0 *
           (1.0 + (0.5 * inv_lambda) + (19.0 / 30.0 * inv_lambda * inv_lambda)));
      return stirling_val + (lambda * (std::log(lambda) - 1.0));
    }

    // Exact series: iterate the Poisson-weighted sum
    //   sum_n accumulates log(n!) = log(2) + log(3) + ... + log(n)
    //   scaled  tracks λ^n / n! (unnormalized Poisson PMF)
    //   accum  = Σ scaled × sum_n
    //   Convergence: stop when new term < 1e-9 × running total
    f64 accum = 0.0;
    f64 sum_n = 0.0;
    f64 scaled = lambda;
    for (int n = 2; n <= 10'000; ++n) {
      sum_n += std::log(static_cast<f64>(n));  // sum_n = log(n!)
      scaled *= lambda / n;                    // scaled = λ^n / n!
      f64 const zscore = scaled * sum_n;
      if (zscore < accum * 1e-9) {
        break;  // convergence check when contribution is negligible
      }
      accum += zscore;
    }
    return accum * std::exp(-lambda);
  }

  // ============================================================================
  // ComputeF: total expected log-factorial for a sequence with ℓ k-mers
  //
  // GC-bias corrected version:
  //   f(ℓ, g) = Σ_{c=0}^k [ C(k,c) · 2^k ] · f_single(ℓ · q_c)
  //   where q_c = (g/2)^c · ((1-g)/2)^{k-c}
  //
  // Groups 4^k possible k-mers into k+1 binomial equivalence classes based
  // on GC content. Each class has C(k,c)·2^k distinct k-mers (c positions
  // chosen for G/C, each can be G or C, remaining k-c can be A or T).
  //
  // When g = 0.5 (uniform), this reduces exactly to the original formula:
  //   f(ℓ) = 4^k × f_single(ℓ / 4^k)
  //
  // This is computed during precomputation only — zero runtime overhead.
  //
  // Example (k=7, ℓ=100, g=0.5):
  //   λ = 100/16384 ≈ 0.006  →  f(100) ≈ 0
  //
  // Example (k=7, ℓ=100, g=0.41):
  //   AT-rich k-mers get higher λ (more expected), GC-rich k-mers get
  //   lower λ (less expected). Net effect: the null model accounts for
  //   compositional bias, reducing false-positive scores in AT-rich regions.
  // ============================================================================
  [[nodiscard]] auto ComputeF(int ell) const -> f64 {
    // Fast path: if perfectly uniform (g ≈ 0.5), use the original 1-bin formula
    if (std::abs(mGc - 0.5) < 1e-6) {
      f64 const lambda = static_cast<f64>(ell) / mNumKmers;
      return static_cast<f64>(mNumKmers) * ComputeFSingle(lambda);
    }

    // GC-bias corrected: group k-mers by GC content
    f64 const safe_gc = std::clamp(mGc, 1e-6, 1.0 - 1e-6);
    f64 const p_gc = safe_gc / 2.0;          // probability of G or C
    f64 const p_at = (1.0 - safe_gc) / 2.0;  // probability of A or T
    f64 const two_pow_k = static_cast<f64>(1ULL << mK);
    f64 total_f = 0.0;

    for (int gc_count = 0; gc_count <= mK; ++gc_count) {
      // Binomial coefficient: C(k, gc_count)
      f64 comb = 1.0;
      for (int j = 1; j <= gc_count; ++j) {
        comb *= static_cast<f64>(mK - j + 1) / static_cast<f64>(j);
      }

      // Number of distinct k-mers with exactly gc_count G/C bases: C(k,gc_count) · 2^k
      f64 const num_kmers = comb * two_pow_k;

      // Probability of one specific k-mer from this composition group
      f64 const prob_kmer = std::pow(p_gc, gc_count) * std::pow(p_at, mK - gc_count);

      // Expected count (λ) for this k-mer
      f64 const lambda = static_cast<f64>(ell) * prob_kmer;

      total_f += num_kmers * ComputeFSingle(lambda);
    }
    return total_f;
  }

  /// Precompute f(ℓ) for all lengths 0..max_len to avoid recomputation.
  void PrecomputeF(int max_len) {
    mF.resize(max_len + 1, 0.0);
    for (int ell = 1; ell <= max_len; ++ell) {
      mF[ell] = ComputeF(ell);
    }
  }
};

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_LONGDUST_SCORER_H_
