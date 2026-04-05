#include "lancet/base/sequence_complexity.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include "absl/strings/str_cat.h"

namespace lancet::base {

// ============================================================================
// SequenceComplexityScorer — Constructor
// ============================================================================

SequenceComplexityScorer::SequenceComplexityScorer(f64 gc_frac)
    : mFlankScorer(FLANK_K, /*max_len=*/1024, gc_frac),
      mHaplotypeScorer(HAPLOTYPE_K, /*max_len=*/4096, gc_frac) {}

// ============================================================================
// ExtractFlank — clamped flanking substring extraction
// ============================================================================

auto SequenceComplexityScorer::ExtractFlank(std::string_view haplotype, const usize var_pos,
                                            const usize var_len, const i64 flank_size) -> std::string_view {
  const auto hap_len = static_cast<i64>(haplotype.size());
  const auto start = std::max(static_cast<i64>(0), static_cast<i64>(var_pos) - flank_size);
  const auto end = std::min(hap_len, static_cast<i64>(var_pos + var_len) + flank_size);
  if (start >= end) return {};
  return haplotype.substr(static_cast<usize>(start), static_cast<usize>(end - start));
}

// ============================================================================
// MaxHomopolymerRun — longest run of identical bases
// ============================================================================

auto SequenceComplexityScorer::MaxHomopolymerRun(std::string_view seq) -> i32 {
  if (seq.empty()) return 0;
  i32 max_run = 1;
  i32 current_run = 1;
  for (usize i = 1; i < seq.size(); ++i) {
    if (seq[i] == seq[i - 1]) {
      current_run++;
      max_run = std::max(max_run, current_run);
    } else {
      current_run = 1;
    }
  }
  return max_run;
}

// ============================================================================
// LocalShannonEntropy — base composition entropy (bits, 0.0–2.0)
//
// H = -Σ p_i log2(p_i) for i ∈ {A, C, G, T}
//
// Interpretation:
//   0.0  = single base type (homopolymer)
//   1.0  = two equally frequent bases
//   2.0  = perfectly balanced ACGT (maximum complexity)
// ============================================================================

auto SequenceComplexityScorer::LocalShannonEntropy(std::string_view seq) -> f32 {
  if (seq.empty()) return 0.0f;

  std::array<usize, 4> counts = {};
  for (const char ch : seq) {
    switch (ch) {
      case 'A': case 'a': counts[0]++; break;
      case 'C': case 'c': counts[1]++; break;
      case 'G': case 'g': counts[2]++; break;
      case 'T': case 't': counts[3]++; break;
      default: break;  // N or other → ignored
    }
  }

  const auto total = static_cast<f32>(counts[0] + counts[1] + counts[2] + counts[3]);
  if (total <= 0.0f) return 0.0f;

  f32 entropy = 0.0f;
  for (const usize cnt : counts) {
    if (cnt == 0) continue;
    const f32 freq = static_cast<f32>(cnt) / total;
    entropy -= freq * std::log2(freq);
  }
  return entropy;
}

// ============================================================================
// IsPrimitiveMotif — check if a motif is a simple period reduction
//
// Example: "ATAT" is NOT primitive (reduces to "AT")
//          "AT" IS primitive (no shorter period)
//          "ABC" IS primitive
// ============================================================================

auto SequenceComplexityScorer::IsPrimitiveMotif(std::string_view motif) -> bool {
  const auto len = static_cast<i32>(motif.size());
  for (i32 p = 1; p < len; ++p) {
    if (len % p != 0) continue;
    bool all_match = true;
    for (i32 i = p; i < len; ++i) {
      if (motif[i] != motif[i % p]) {
        all_match = false;
        break;
      }
    }
    if (all_match) return false;  // reducible to period p
  }
  return true;
}

// ============================================================================
// FindExactRepeats — scan for exact tandem repeats
//
// For each position and each period 1..max_period, check if the motif
// starting at that position repeats ≥ min_copies times.
//
// Filters non-primitive motifs (e.g., ATAT when AT is found).
// ============================================================================

auto SequenceComplexityScorer::FindExactRepeats(std::string_view seq,
                                                const i32 max_period,
                                                const f32 min_copies) -> std::vector<TandemRepeatResult> {
  std::vector<TandemRepeatResult> results;
  const auto seq_len = static_cast<i32>(seq.size());

  for (i32 period = 1; period <= max_period && period <= seq_len; ++period) {
    for (i32 start = 0; start <= seq_len - period; ++start) {
      const auto motif = seq.substr(start, period);

      // Skip non-primitive motifs (e.g., ATAT when AT covers it)
      if (period > 1 && !IsPrimitiveMotif(motif)) continue;

      // Count consecutive exact matches
      i32 match_len = period;  // first copy
      while (start + match_len + period <= seq_len) {
        bool matches = true;
        for (i32 j = 0; j < period; ++j) {
          if (seq[start + match_len + j] != motif[j]) {
            matches = false;
            break;
          }
        }
        if (!matches) break;
        match_len += period;
      }

      // Check partial copy at the end
      i32 partial = 0;
      while (start + match_len + partial < seq_len && partial < period &&
             seq[start + match_len + partial] == motif[partial]) {
        partial++;
      }

      const f32 copies = static_cast<f32>(match_len + partial) / static_cast<f32>(period);
      if (copies >= min_copies) {
        results.push_back({
            .period = period,
            .copies = copies,
            .start_pos = start,
            .span_length = match_len + partial,
            .total_errors = 0,
            .is_exact = true,
        });
        // Skip ahead past this repeat to avoid overlapping reports
        start += match_len - 1;
      }
    }
  }
  return results;
}

// ============================================================================
// FindApproxRepeats — scan for approximate tandem repeats
//
// Uses banded edit-distance matching: for each position and period,
// check if successive copies match the motif within max_edits_per_unit.
// Only reports results with purity ≥ 0.75.
// ============================================================================

auto SequenceComplexityScorer::FindApproxRepeats(std::string_view seq,
                                                 const i32 max_period,
                                                 const f32 min_copies,
                                                 const i32 max_edits_per_unit) -> std::vector<TandemRepeatResult> {
  std::vector<TandemRepeatResult> results;
  const auto seq_len = static_cast<i32>(seq.size());

  for (i32 period = 1; period <= max_period && period <= seq_len; ++period) {
    for (i32 start = 0; start <= seq_len - period; ++start) {
      const auto motif = seq.substr(start, period);
      if (period > 1 && !IsPrimitiveMotif(motif)) continue;

      i32 total_span = period;
      i32 total_errors = 0;
      i32 num_units = 1;

      // Extend by comparing successive units to the motif
      while (start + total_span + period <= seq_len) {
        i32 unit_errors = 0;
        for (i32 j = 0; j < period; ++j) {
          if (seq[start + total_span + j] != motif[j]) {
            unit_errors++;
          }
        }
        if (unit_errors > max_edits_per_unit) break;
        total_errors += unit_errors;
        total_span += period;
        num_units++;
      }

      const f32 copies = static_cast<f32>(total_span) / static_cast<f32>(period);
      const f32 purity = total_span > 0
          ? 1.0f - static_cast<f32>(total_errors) / static_cast<f32>(total_span)
          : 0.0f;

      if (copies >= min_copies && purity >= 0.75f) {
        results.push_back({
            .period = period,
            .copies = copies,
            .start_pos = start,
            .span_length = total_span,
            .total_errors = total_errors,
            .is_exact = false,
        });
        start += total_span - 1;
      }
    }
  }
  return results;
}

// ============================================================================
// FlattenTRFeatures — collapse TR hits into fixed-width features
//
// Retains only the fields consumed by the SequenceComplexity distiller:
//   dist_to_nearest_tr, nearest_tr_period, nearest_tr_purity, is_stutter_indel
// ============================================================================

auto SequenceComplexityScorer::FlattenTRFeatures(
    const std::vector<TandemRepeatResult>& results,
    const i32 variant_pos, const i32 variant_length,
    const i32 /*window_size*/) -> VariantTRFeatures {

  VariantTRFeatures feat;
  if (results.empty()) return feat;

  i32 nearest_dist = std::numeric_limits<i32>::max();
  const i32 var_end = variant_pos + variant_length;

  for (const auto& tr : results) {
    const i32 tr_end = tr.start_pos + tr.span_length;

    // Distance from variant to this TR
    i32 dist = 0;
    if (variant_pos >= tr.start_pos && variant_pos < tr_end) {
      dist = 0;  // inside
    } else if (variant_pos < tr.start_pos) {
      dist = tr.start_pos - var_end;
    } else {
      dist = variant_pos - tr_end;
    }
    dist = std::max(0, dist);

    if (dist < nearest_dist) {
      nearest_dist = dist;
      feat.dist_to_nearest_tr = dist;
      feat.nearest_tr_period = tr.period;
      feat.nearest_tr_purity = tr.Purity();
    }

    // Check for stutter indel pattern
    if (dist <= 1 && variant_length > 0 && variant_length <= tr.period) {
      feat.is_stutter_indel = 1;
    }
  }

  return feat;
}

// ============================================================================
// AccumulateTRFeatures — run motif detection and take element-wise max
// ============================================================================

void SequenceComplexityScorer::AccumulateTRFeatures(VariantTRFeatures& features,
                                                    std::string_view window,
                                                    const i32 var_pos_in_window,
                                                    const i32 var_length,
                                                    const i32 window_size) {
  auto exact = FindExactRepeats(window);
  auto approx = FindApproxRepeats(window);

  // Merge results
  std::vector<TandemRepeatResult> all_results;
  all_results.reserve(exact.size() + approx.size());
  all_results.insert(all_results.end(), exact.begin(), exact.end());
  all_results.insert(all_results.end(), approx.begin(), approx.end());

  auto new_feat = FlattenTRFeatures(all_results, var_pos_in_window, var_length, window_size);

  // Take the closest TR (minimum distance)
  if (new_feat.dist_to_nearest_tr >= 0 &&
      (features.dist_to_nearest_tr < 0 || new_feat.dist_to_nearest_tr < features.dist_to_nearest_tr)) {
    features.dist_to_nearest_tr = new_feat.dist_to_nearest_tr;
    features.nearest_tr_period = new_feat.nearest_tr_period;
    features.nearest_tr_purity = new_feat.nearest_tr_purity;
  }
  features.is_stutter_indel = std::max(features.is_stutter_indel, new_feat.is_stutter_indel);
}

// ============================================================================
// Score — main entry point. Delegates to ScoreContext, ScoreDeltas, ScoreTrMotif.
// ============================================================================

auto SequenceComplexityScorer::Score(
    std::string_view ref_haplotype, const usize ref_pos, const usize ref_len,
    std::string_view alt_haplotype, const usize alt_pos, const usize alt_len) const
    -> SequenceComplexity {
  SequenceComplexity cx;
  ScoreContext(cx, ref_haplotype, ref_pos, ref_len);
  ScoreDeltas(cx, ref_haplotype, ref_pos, ref_len, alt_haplotype, alt_pos, alt_len);
  ScoreTrMotif(cx, alt_haplotype, alt_pos, alt_len);
  return cx;
}

// ============================================================================
// ScoreContext — "How brittle is the genome here?" (strictly REF)
//
// Populates: mContextHRun, mContextEntropy, mContextFlankLQ, mContextHaplotypeLQ
// ============================================================================

void SequenceComplexityScorer::ScoreContext(SequenceComplexity& cx,
                                            std::string_view ref_hap,
                                            const usize ref_pos,
                                            const usize ref_len) const {
  // HRun + Entropy at ±20bp
  const auto ctx_window = ExtractFlank(ref_hap, ref_pos, ref_len, CONTEXT_FLANK);
  cx.mContextHRun = MaxHomopolymerRun(ctx_window);
  cx.mContextEntropy = LocalShannonEntropy(ctx_window);

  // LongdustQ (k=4) at ±50bp — log1p-squashed to compress heavy tails
  const auto lq_window = ExtractFlank(ref_hap, ref_pos, ref_len, LQ_FLANK);
  cx.mContextFlankLQ = std::log1p(std::max(0.0, mFlankScorer.Score(lq_window)));

  // LongdustQ (k=7) on full haplotype — log1p-squashed
  cx.mContextHaplotypeLQ = std::log1p(std::max(0.0, mHaplotypeScorer.Score(ref_hap)));
}

// ============================================================================
// ScoreDeltas — "How did the variant alter complexity?" (ALT minus REF)
//
// Populates: mDeltaHRun, mDeltaEntropy, mDeltaFlankLQ
//
// Context + Perturbation paradigm: REF context features tell the model
// "the genome is brittle here", while deltas tell it "the variant made it
// worse/better". An EBM can learn rescue logic: even if context is dangerous,
// a negative delta (variant broke the homopolymer) rescues the call.
// ============================================================================

void SequenceComplexityScorer::ScoreDeltas(SequenceComplexity& cx,
                                           std::string_view ref_hap, const usize ref_pos, const usize ref_len,
                                           std::string_view alt_hap, const usize alt_pos, const usize alt_len) const {
  // HRun delta at ±5bp
  const auto ref_hrun_window = ExtractFlank(ref_hap, ref_pos, ref_len, DELTA_HRUN_FLANK);
  const auto alt_hrun_window = ExtractFlank(alt_hap, alt_pos, alt_len, DELTA_HRUN_FLANK);
  cx.mDeltaHRun = MaxHomopolymerRun(alt_hrun_window) - MaxHomopolymerRun(ref_hrun_window);

  // Entropy delta at ±10bp
  const auto ref_ent_window = ExtractFlank(ref_hap, ref_pos, ref_len, DELTA_ENTROPY_FLANK);
  const auto alt_ent_window = ExtractFlank(alt_hap, alt_pos, alt_len, DELTA_ENTROPY_FLANK);
  cx.mDeltaEntropy = LocalShannonEntropy(alt_ent_window) - LocalShannonEntropy(ref_ent_window);

  // LongdustQ delta at ±50bp (log-space)
  const auto alt_lq_window = ExtractFlank(alt_hap, alt_pos, alt_len, LQ_FLANK);
  const auto alt_lq = std::log1p(std::max(0.0, mFlankScorer.Score(alt_lq_window)));
  cx.mDeltaFlankLQ = alt_lq - cx.mContextFlankLQ;
}

// ============================================================================
// ScoreTrMotif — "What is the repeat environment of the ALT allele?"
//
// Populates: mTrAffinity, mTrPurity, mTrPeriod, mIsStutterIndel
//
// All features are sentinel-safe: when no TR is found (dist_to_nearest_tr < 0),
// the affinity/purity/period are set to 0 rather than propagating the -1
// sentinel, which would break EBM numerical binning.
// ============================================================================

void SequenceComplexityScorer::ScoreTrMotif(SequenceComplexity& cx,
                                            std::string_view alt_hap,
                                            const usize alt_pos,
                                            const usize alt_len) const {
  const auto window = ExtractFlank(alt_hap, alt_pos, alt_len, TR_MOTIF_FLANK);
  const auto start = std::max(static_cast<i64>(0), static_cast<i64>(alt_pos) - TR_MOTIF_FLANK);
  const auto var_pos_in_window = static_cast<i32>(static_cast<i64>(alt_pos) - start);

  VariantTRFeatures tr;
  AccumulateTRFeatures(tr, window, var_pos_in_window,
                       static_cast<i32>(alt_len), static_cast<i32>(window.size()));

  // Sentinel-safe transforms: dist < 0 means no TR found → all zeros
  if (tr.dist_to_nearest_tr < 0) {
    cx.mTrAffinity = 0.0f;
    cx.mTrPurity = 0.0f;
    cx.mTrPeriod = 0;
  } else {
    cx.mTrAffinity = 1.0f / (1.0f + static_cast<f32>(tr.dist_to_nearest_tr));
    cx.mTrPurity = tr.nearest_tr_purity;
    cx.mTrPeriod = tr.nearest_tr_period;
  }
  cx.mIsStutterIndel = tr.is_stutter_indel;
}

// ============================================================================
// SequenceComplexity::FormatVcfValue — 11 comma-separated values
//
// Order: Context(4), Delta(3), TrMotif(4)
// ============================================================================

auto SequenceComplexity::FormatVcfValue() const -> std::string {
  return absl::StrCat(
      mContextHRun, ",",
      FormatComplexityScore(mContextEntropy), ",",
      FormatComplexityScore(mContextFlankLQ), ",",
      FormatComplexityScore(mContextHaplotypeLQ), ",",
      mDeltaHRun, ",",
      FormatComplexityScore(mDeltaEntropy), ",",
      FormatComplexityScore(mDeltaFlankLQ), ",",
      FormatComplexityScore(mTrAffinity), ",",
      FormatComplexityScore(mTrPurity), ",",
      mTrPeriod, ",",
      mIsStutterIndel);
}

// ============================================================================
// SequenceComplexity::MergeMax — element-wise max for multi-haplotype
// and multi-allelic merging.
//
// Context features: max gives highest brittleness.
// Delta features: max gives most extreme perturbation.
// TR motif features: max gives closest/purest repeat.
// ============================================================================

void SequenceComplexity::MergeMax(const SequenceComplexity& other) {
  // Context
  mContextHRun = std::max(mContextHRun, other.mContextHRun);
  mContextEntropy = std::max(mContextEntropy, other.mContextEntropy);
  mContextFlankLQ = std::max(mContextFlankLQ, other.mContextFlankLQ);
  mContextHaplotypeLQ = std::max(mContextHaplotypeLQ, other.mContextHaplotypeLQ);

  // Deltas
  mDeltaHRun = std::max(mDeltaHRun, other.mDeltaHRun);
  mDeltaEntropy = std::max(mDeltaEntropy, other.mDeltaEntropy);
  mDeltaFlankLQ = std::max(mDeltaFlankLQ, other.mDeltaFlankLQ);

  // TR Motif
  mTrAffinity = std::max(mTrAffinity, other.mTrAffinity);
  mTrPurity = std::max(mTrPurity, other.mTrPurity);
  mTrPeriod = std::max(mTrPeriod, other.mTrPeriod);
  mIsStutterIndel = std::max(mIsStutterIndel, other.mIsStutterIndel);
}

}  // namespace lancet::base
