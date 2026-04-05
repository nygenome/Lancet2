#ifndef SRC_LANCET_BASE_SEQUENCE_COMPLEXITY_H_
#define SRC_LANCET_BASE_SEQUENCE_COMPLEXITY_H_

#include <cmath>
#include <string>
#include <string_view>
#include <vector>

#include "absl/strings/str_cat.h"
#include "lancet/base/longdust_scorer.h"
#include "lancet/base/types.h"

namespace lancet::base {

// ============================================================================
// Tandem Repeat Detection — Result and Intermediate Types
// ============================================================================

/// Result of detecting a single tandem repeat within a sequence window.
struct TandemRepeatResult {
  i32 period = 0;           ///< motif length (1=homopolymer, 2=dinucleotide, ...)
  f32 copies = 0.0f;        ///< fractional copy count (span / period)
  i32 start_pos = 0;        ///< 0-based start within scored window
  i32 span_length = 0;      ///< total bp covered by the repeat
  i32 total_errors = 0;     ///< edit distance sum (approximate only, 0 for exact)
  bool is_exact = false;    ///< true if found by exact matching

  /// Fraction of span that is error-free.
  /// 1.0 for exact repeats. For approximate: 1.0 - total_errors / span_length.
  [[nodiscard]] auto Purity() const -> f32 {
    if (span_length <= 0) return 0.0f;
    return 1.0f - static_cast<f32>(total_errors) / static_cast<f32>(span_length);
  }
};

/// Intermediate tandem repeat features for a single flanking window.
/// Used internally by SequenceComplexityScorer during motif detection.
/// Only the fields consumed by the SequenceComplexity distiller are retained.
struct VariantTRFeatures {
  i32 dist_to_nearest_tr = -1;  ///< bp distance to nearest TR (-1 = no TR found)
  i32 nearest_tr_period = 0;    ///< period of nearest TR
  f32 nearest_tr_purity = 0.0f; ///< purity of nearest TR
  i32 is_stutter_indel = 0;     ///< 1 if INDEL length ≤ period and adjacent to TR
};

// ============================================================================
// SequenceComplexity — 11-feature ML-ready sequence complexity vector
//
// Distills raw multi-scale sequence metrics (homopolymer runs, Shannon
// entropy, LongdustQ k-mer concentration, tandem repeat motifs) into an
// orthogonal feature set designed for additive ML models (EBMs/GAMs).
//
// Three conceptual groups:
//
// ── Context (4 features, strictly REF) ──────────────────────────────
//   "How brittle is the genome here, regardless of the variant?"
//
//   ContextHRun:        Max homopolymer run in REF ±20bp.
//   ContextEntropy:     Shannon entropy (bits) in REF ±20bp.
//   ContextFlankLQ:     log1p-squashed LongdustQ (k=4) in REF ±50bp.
//   ContextHaplotypeLQ: log1p-squashed LongdustQ (k=7) on full REF haplotype.
//
// ── Deltas (3 features, ALT minus REF) ──────────────────────────────
//   "How did the variant alter the local sequence complexity?"
//
//   DeltaHRun:    ALT ±5bp HRun − REF ±5bp HRun.
//                 Positive = variant extended a homopolymer (artifact signal).
//                 Negative = variant broke a homopolymer (rescue signal).
//   DeltaEntropy: ALT ±10bp entropy − REF ±10bp entropy.
//                 Negative = variant reduced diversity (gap mimicking deletion).
//   DeltaFlankLQ: log1p(ALT ±50bp LQ) − log1p(REF ±50bp LQ).
//                 Positive = variant exacerbated microsatellite repetitiveness.
//
// ── TR Motif (4 features, strictly ALT ±50bp) ───────────────────────
//   "What is the tandem repeat environment of the final allele?"
//
//   TrAffinity:     1/(1+dist) where dist = distance to nearest TR.
//                   Sentinel-safe: dist<0 → 0 (no TR found), dist=0 → 1 (on TR).
//   TrPurity:       Purity of nearest TR. dist<0 → 0.
//   TrPeriod:       Period of nearest TR. dist<0 → 0.
//   IsStutterIndel: 1 if INDEL length ≤ period and adjacent to TR, 0 otherwise.
//
// ============================================================================
class SequenceComplexity {
 public:
  // ── Context accessors ──────────────────────────────────────────────
  [[nodiscard]] auto ContextHRun() const noexcept -> i32 { return mContextHRun; }
  [[nodiscard]] auto ContextEntropy() const noexcept -> f32 { return mContextEntropy; }
  [[nodiscard]] auto ContextFlankLQ() const noexcept -> f64 { return mContextFlankLQ; }
  [[nodiscard]] auto ContextHaplotypeLQ() const noexcept -> f64 { return mContextHaplotypeLQ; }

  // ── Delta accessors ────────────────────────────────────────────────
  [[nodiscard]] auto DeltaHRun() const noexcept -> i32 { return mDeltaHRun; }
  [[nodiscard]] auto DeltaEntropy() const noexcept -> f32 { return mDeltaEntropy; }
  [[nodiscard]] auto DeltaFlankLQ() const noexcept -> f64 { return mDeltaFlankLQ; }

  // ── TR Motif accessors ─────────────────────────────────────────────
  [[nodiscard]] auto TrAffinity() const noexcept -> f32 { return mTrAffinity; }
  [[nodiscard]] auto TrPurity() const noexcept -> f32 { return mTrPurity; }
  [[nodiscard]] auto TrPeriod() const noexcept -> i32 { return mTrPeriod; }
  [[nodiscard]] auto IsStutterIndel() const noexcept -> i32 { return mIsStutterIndel; }

  /// Format as 11 comma-separated values for VCF SEQ_CX INFO tag.
  [[nodiscard]] auto FormatVcfValue() const -> std::string;

  /// Element-wise max across two SequenceComplexity objects.
  /// Used for multi-haplotype and multi-allelic merging.
  void MergeMax(const SequenceComplexity& other);

 private:
  // Memory alignment: 8B → 4B → 4B (descending alignment convention)
  f64 mContextFlankLQ = 0.0;      ///< log1p-squashed REF ±50bp LQ (k=4)
  f64 mContextHaplotypeLQ = 0.0;  ///< log1p-squashed REF full-haplotype LQ (k=7)
  f64 mDeltaFlankLQ = 0.0;        ///< log-space ALT−REF delta at ±50bp
  f32 mContextEntropy = 0.0f;     ///< Shannon entropy in REF ±20bp
  f32 mDeltaEntropy = 0.0f;       ///< ALT−REF entropy delta at ±10bp
  f32 mTrAffinity = 0.0f;         ///< sentinel-safe [0,1]
  f32 mTrPurity = 0.0f;           ///< sentinel-safe [0,1]
  i32 mContextHRun = 0;           ///< max homopolymer run in REF ±20bp
  i32 mDeltaHRun = 0;             ///< ALT−REF HRun delta at ±5bp
  i32 mTrPeriod = 0;              ///< sentinel-safe (0 = no TR)
  i32 mIsStutterIndel = 0;        ///< binary stutter flag

  friend class SequenceComplexityScorer;  // populates private fields
};

// ============================================================================
// SequenceComplexityScorer — produces SequenceComplexity from haplotypes
//
// Bundles all sequence complexity methods (homopolymer run, Shannon entropy,
// tandem repeat motif detection, LongdustQ scoring) and produces the
// distilled 11-feature SequenceComplexity output via Score().
//
// Internally organized by feature group (context / delta / tr_motif),
// not by scale (ultra / micro / macro). Each group extracts the specific
// flanking window it needs and computes only the features it contributes.
//
// The class owns its LongdustQScorer instances internally and is constructed
// once in VariantBuilder, shared across all windows (same thread lifecycle).
//
// Motif Detection Configuration:
//   max_period = 6: covers homopolymers through hexanucleotide repeats.
//     Period 1–3 repeats account for >90% of Illumina INDEL errors.
//
//   min_copies_exact = 2.5: requires ≥2.5 full copies for exact matches.
//     Below this, a "repeat" is just two overlapping instances.
//
//   min_copies_approx = 3.0: higher threshold for approximate matches
//     to compensate for relaxed matching.
//
//   max_edits_per_unit = 1: allows 1 mismatch/indel per repeat unit.
//
//   min_purity = 0.75: fraction of span that is error-free. Below 0.75,
//     the "repeat" is too degraded to cause polymerase stutter.
// ============================================================================
class SequenceComplexityScorer {
 public:
  /// @param gc_frac  Global background GC fraction for LongdustQ scoring.
  ///                 Default: 0.41 (human genome-wide average).
  ///                 Set to 0.5 for uniform distribution (no GC correction).
  ///                 See --genome-gc-bias CLI parameter.
  explicit SequenceComplexityScorer(f64 gc_frac = 0.41);

  /// Score sequence complexity for a single variant.
  /// Produces the 11-feature SequenceComplexity by delegating to
  /// ScoreContext, ScoreDeltas, and ScoreTrMotif.
  ///
  /// @param ref_haplotype  Full assembled reference haplotype sequence
  /// @param ref_pos        0-based variant start position in REF haplotype
  /// @param ref_len        REF allele length
  /// @param alt_haplotype  Full assembled ALT haplotype sequence
  /// @param alt_pos        0-based variant start position in ALT haplotype
  /// @param alt_len        Max of REF/ALT allele lengths
  [[nodiscard]] auto Score(
      std::string_view ref_haplotype, usize ref_pos, usize ref_len,
      std::string_view alt_haplotype, usize alt_pos, usize alt_len) const
      -> SequenceComplexity;

  // ── Component methods (public for unit testing) ───────────────────────

  /// Maximum homopolymer run length in a sequence.
  /// Returns 0 for empty sequences.
  [[nodiscard]] static auto MaxHomopolymerRun(std::string_view seq) -> i32;

  /// Shannon entropy of base frequencies in a sequence (bits, 0.0–2.0).
  /// H = -Σ p_i log2(p_i) where p_i is frequency of base i (A/C/G/T).
  /// Returns 0.0 for sequences with only one base type.
  [[nodiscard]] static auto LocalShannonEntropy(std::string_view seq) -> f32;

  // ── Motif detection ───────────────────────────────────────────────────

  /// Find exact tandem repeats by checking if any motif of period 1..max_period
  /// repeats ≥ min_copies times starting at each position.
  /// Filters primitive motifs (e.g., ATAT collapsed to AT).
  [[nodiscard]] static auto FindExactRepeats(std::string_view seq,
                                             i32 max_period = 6,
                                             f32 min_copies = 2.5f)
      -> std::vector<TandemRepeatResult>;

  /// Find approximate tandem repeats using edit-distance matching.
  /// Each repeat unit may have up to max_edits_per_unit errors.
  /// Only reports results with purity ≥ 0.75.
  [[nodiscard]] static auto FindApproxRepeats(std::string_view seq,
                                              i32 max_period = 6,
                                              f32 min_copies = 3.0f,
                                              i32 max_edits_per_unit = 1)
      -> std::vector<TandemRepeatResult>;

  // ── Feature flattening ────────────────────────────────────────────────

  /// Flatten a list of TR hits into VariantTRFeatures relative to a variant
  /// at (variant_pos, variant_length) within a window of total size window_size.
  [[nodiscard]] static auto FlattenTRFeatures(
      const std::vector<TandemRepeatResult>& results,
      i32 variant_pos, i32 variant_length, i32 window_size) -> VariantTRFeatures;

 private:
  // Owned LongdustQ scorers — constructed once, used across all windows.
  LongdustQScorer mFlankScorer;         ///< k=4 for ±50bp flanks
  LongdustQScorer mHaplotypeScorer;     ///< k=7 for full haplotype

  // LongdustQ k-mer sizes
  static constexpr int FLANK_K = 4;
  static constexpr int HAPLOTYPE_K = 7;

  // ── Flank sizes organized by feature group ────────────────────────────
  static constexpr i64 CONTEXT_FLANK = 20;       ///< HRun + Entropy context
  static constexpr i64 DELTA_HRUN_FLANK = 5;      ///< HRun delta
  static constexpr i64 DELTA_ENTROPY_FLANK = 10;  ///< Entropy delta
  static constexpr i64 LQ_FLANK = 50;             ///< LongdustQ (context + delta)
  static constexpr i64 TR_MOTIF_FLANK = 50;        ///< TR motif detection

  // ── Private scoring helpers (one per feature group) ───────────────────

  /// Populate context features from REF haplotype.
  void ScoreContext(SequenceComplexity& cx,
                    std::string_view ref_hap, usize ref_pos, usize ref_len) const;

  /// Populate delta features from ALT−REF differences.
  void ScoreDeltas(SequenceComplexity& cx,
                   std::string_view ref_hap, usize ref_pos, usize ref_len,
                   std::string_view alt_hap, usize alt_pos, usize alt_len) const;

  /// Populate TR motif features from ALT haplotype.
  void ScoreTrMotif(SequenceComplexity& cx,
                    std::string_view alt_hap, usize alt_pos, usize alt_len) const;

  // ── Internal helpers ──────────────────────────────────────────────────

  /// Returns true if the motif is a simple repeat of a shorter primitive.
  [[nodiscard]] static auto IsPrimitiveMotif(std::string_view motif) -> bool;

  /// Extract a clamped flanking substring from haplotype around var_pos.
  [[nodiscard]] static auto ExtractFlank(std::string_view haplotype, usize var_pos,
                                         usize var_len, i64 flank_size) -> std::string_view;

  /// Run motif detection (both exact + approx) on a flanking window and
  /// take element-wise max into existing VariantTRFeatures.
  static void AccumulateTRFeatures(VariantTRFeatures& features,
                                   std::string_view window,
                                   i32 var_pos_in_window,
                                   i32 var_length,
                                   i32 window_size);
};

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_SEQUENCE_COMPLEXITY_H_
