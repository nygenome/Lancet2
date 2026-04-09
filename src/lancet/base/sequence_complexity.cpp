#include "lancet/base/sequence_complexity.h"

#include "lancet/base/longdust_scorer.h"
#include "lancet/base/types.h"

#include "absl/strings/str_cat.h"

#include <algorithm>
#include <array>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include <cmath>

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

auto SequenceComplexityScorer::ExtractFlank(std::string_view haplotype, usize const var_pos,
                                            usize const var_len, i64 const flank_size)
    -> std::string_view {
  auto const hap_len = static_cast<i64>(haplotype.size());
  auto const start = std::max(static_cast<i64>(0), static_cast<i64>(var_pos) - flank_size);
  auto const end = std::min(hap_len, static_cast<i64>(var_pos + var_len) + flank_size);
  if (start >= end) {
    return {};
  }
  return haplotype.substr(static_cast<usize>(start), static_cast<usize>(end - start));
}

// ============================================================================
// MaxHomopolymerRun — longest run of identical bases
// ============================================================================

auto SequenceComplexityScorer::MaxHomopolymerRun(std::string_view seq) -> i32 {
  if (seq.empty()) {
    return 0;
  }
  i32 max_run = 1;
  i32 current_run = 1;
  for (usize i = 1; i < seq.size(); ++i) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
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
  if (seq.empty()) {
    return 0.0F;
  }

  std::array<usize, 4> counts = {};
  for (char const chr : seq) {
    switch (chr) {
      // NOLINTBEGIN(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
      case 'A':
      case 'a':
        counts[0]++;
        break;
      case 'C':
      case 'c':
        counts[1]++;
        break;
      case 'G':
      case 'g':
        counts[2]++;
        break;
      case 'T':
      case 't':
        counts[3]++;
        break;
      // NOLINTEND(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
      default:
        break;  // N or other → ignored
    }
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
  auto const total = static_cast<f32>(counts[0] + counts[1] + counts[2] + counts[3]);
  if (total <= 0.0F) {
    return 0.0F;
  }

  f32 entropy = 0.0F;
  for (usize const cnt : counts) {
    if (cnt == 0) {
      continue;
    }
    f32 const freq = static_cast<f32>(cnt) / total;
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
  auto const len = static_cast<i32>(motif.size());
  for (i32 period = 1; period < len; ++period) {
    if (len % period != 0) {
      continue;
    }
    bool all_match = true;
    for (i32 i = period; i < len; ++i) {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
      if (motif[i] != motif[i % period]) {
        all_match = false;
        break;
      }
    }
    if (all_match) {
      return false;  // reducible to period p
    }
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

// NOLINTNEXTLINE(readability-function-size)  // TODO(lancet): refactor to reduce function size
auto SequenceComplexityScorer::FindExactRepeats(std::string_view seq, i32 const max_period,
                                                f32 const min_copies)
    -> std::vector<TandemRepeatResult> {
  std::vector<TandemRepeatResult> results;
  auto const seq_len = static_cast<i32>(seq.size());

  for (i32 period = 1; period <= max_period && period <= seq_len; ++period) {
    for (i32 start = 0; start <= seq_len - period; ++start) {
      auto const motif = seq.substr(start, period);

      // Skip non-primitive motifs (e.g., ATAT when AT covers it)
      if (period > 1 && !IsPrimitiveMotif(motif)) {
        continue;
      }

      // Count consecutive exact matches
      i32 match_len = period;  // first copy
      while (start + match_len + period <= seq_len) {
        bool matches = true;
        for (i32 j = 0; j < period; ++j) {
          // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
          if (seq[start + match_len + j] != motif[j]) {
            matches = false;
            break;
          }
        }
        if (!matches) {
          break;
        }
        match_len += period;
      }

      // Check partial copy at the end
      i32 partial = 0;
      while (start + match_len + partial < seq_len &&
             partial < period &&
             // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
             seq[start + match_len + partial] == motif[partial]) {
        partial++;
      }

      f32 const copies = static_cast<f32>(match_len + partial) / static_cast<f32>(period);
      if (copies >= min_copies) {
        results.push_back({
            .mPeriod = period,
            .mCopies = copies,
            .mStartPos = start,
            .mSpanLength = match_len + partial,
            .mTotalErrors = 0,
            .mIsExact = true,
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

// NOLINTNEXTLINE(readability-function-size)  // TODO(lancet): refactor to reduce function size
auto SequenceComplexityScorer::FindApproxRepeats(std::string_view seq, i32 const max_period,
                                                 f32 const min_copies, i32 const max_edits_per_unit)
    -> std::vector<TandemRepeatResult> {
  std::vector<TandemRepeatResult> results;
  auto const seq_len = static_cast<i32>(seq.size());

  for (i32 period = 1; period <= max_period && period <= seq_len; ++period) {
    for (i32 start = 0; start <= seq_len - period; ++start) {
      auto const motif = seq.substr(start, period);
      if (period > 1 && !IsPrimitiveMotif(motif)) {
        continue;
      }

      i32 total_span = period;
      i32 total_errors = 0;
      i32 num_units = 1;

      // Extend by comparing successive units to the motif
      while (start + total_span + period <= seq_len) {
        i32 unit_errors = 0;
        for (i32 j = 0; j < period; ++j) {
          // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
          if (seq[start + total_span + j] != motif[j]) {
            unit_errors++;
          }
        }
        if (unit_errors > max_edits_per_unit) {
          break;
        }
        total_errors += unit_errors;
        total_span += period;
        num_units++;
      }

      f32 const copies = static_cast<f32>(total_span) / static_cast<f32>(period);
      f32 const purity =
          total_span > 0 ? 1.0F - (static_cast<f32>(total_errors) / static_cast<f32>(total_span))
                         : 0.0F;

      if (copies >= min_copies && purity >= 0.75F) {
        results.push_back({
            .mPeriod = period,
            .mCopies = copies,
            .mStartPos = start,
            .mSpanLength = total_span,
            .mTotalErrors = total_errors,
            .mIsExact = false,
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

auto SequenceComplexityScorer::FlattenTRFeatures(std::vector<TandemRepeatResult> const& results,
                                                 i32 const variant_pos, i32 const variant_length,
                                                 i32 const /*window_size*/) -> VariantTRFeatures {
  VariantTRFeatures feat;
  if (results.empty()) {
    return feat;
  }

  i32 nearest_dist = std::numeric_limits<i32>::max();
  i32 const var_end = variant_pos + variant_length;

  for (auto const& trep : results) {
    i32 const tr_end = trep.mStartPos + trep.mSpanLength;

    // Distance from variant to this TR
    i32 dist = 0;
    if (variant_pos >= trep.mStartPos && variant_pos < tr_end) {
      dist = 0;  // inside
    } else if (variant_pos < trep.mStartPos) {
      dist = trep.mStartPos - var_end;
    } else {
      dist = variant_pos - tr_end;
    }
    dist = std::max(0, dist);

    if (dist < nearest_dist) {
      nearest_dist = dist;
      feat.mDistToNearestTr = dist;
      feat.mNearestTrPeriod = trep.mPeriod;
      feat.mNearestTrPurity = trep.Purity();
    }

    // Check for stutter indel pattern
    if (dist <= 1 && variant_length > 0 && variant_length <= trep.mPeriod) {
      feat.mIsStutterIndel = 1;
    }
  }

  return feat;
}

// ============================================================================
// AccumulateTRFeatures — run motif detection and take element-wise max
// ============================================================================

void SequenceComplexityScorer::AccumulateTRFeatures(VariantTRFeatures& features,
                                                    std::string_view window,
                                                    i32 const var_pos_in_window,
                                                    i32 const var_length, i32 const window_size) {
  auto exact = FindExactRepeats(window);
  auto approx = FindApproxRepeats(window);

  // Merge results
  std::vector<TandemRepeatResult> all_results;
  all_results.reserve(exact.size() + approx.size());
  all_results.insert(all_results.end(), exact.begin(), exact.end());
  all_results.insert(all_results.end(), approx.begin(), approx.end());

  auto new_feat = FlattenTRFeatures(all_results, var_pos_in_window, var_length, window_size);

  // Take the closest TR (minimum distance)
  if (new_feat.mDistToNearestTr >= 0 &&
      (features.mDistToNearestTr < 0 || new_feat.mDistToNearestTr < features.mDistToNearestTr)) {
    features.mDistToNearestTr = new_feat.mDistToNearestTr;
    features.mNearestTrPeriod = new_feat.mNearestTrPeriod;
    features.mNearestTrPurity = new_feat.mNearestTrPurity;
  }
  features.mIsStutterIndel = std::max(features.mIsStutterIndel, new_feat.mIsStutterIndel);
}

// ============================================================================
// Score — main entry point. Delegates to ScoreContext, ScoreDeltas, ScoreTrMotif.
// ============================================================================

auto SequenceComplexityScorer::Score(std::string_view ref_haplotype, usize const ref_pos,
                                     usize const ref_len, std::string_view alt_haplotype,
                                     usize const alt_pos, usize const alt_len) const
    -> SequenceComplexity {
  SequenceComplexity cplx;
  ScoreContext(cplx, ref_haplotype, ref_pos, ref_len);
  ScoreDeltas(cplx, ref_haplotype, ref_pos, ref_len, alt_haplotype, alt_pos, alt_len);
  ScoreTrMotif(cplx, alt_haplotype, alt_pos, alt_len);
  return cplx;
}

// ============================================================================
// ScoreContext — "How brittle is the genome here?" (strictly REF)
//
// Populates: mContextHRun, mContextEntropy, mContextFlankLQ, mContextHaplotypeLQ
// ============================================================================

void SequenceComplexityScorer::ScoreContext(SequenceComplexity& cplx, std::string_view ref_hap,
                                            usize const ref_pos, usize const ref_len) const {
  // HRun + Entropy at ±20bp
  auto const ctx_window = ExtractFlank(ref_hap, ref_pos, ref_len, CONTEXT_FLANK);
  cplx.mContextHRun = MaxHomopolymerRun(ctx_window);
  cplx.mContextEntropy = LocalShannonEntropy(ctx_window);

  // LongdustQ (k=4) at ±50bp — log1p-squashed to compress heavy tails
  auto const lq_window = ExtractFlank(ref_hap, ref_pos, ref_len, LQ_FLANK);
  cplx.mContextFlankLQ = std::log1p(std::max(0.0, mFlankScorer.Score(lq_window)));

  // LongdustQ (k=7) on full haplotype — log1p-squashed
  cplx.mContextHaplotypeLQ = std::log1p(std::max(0.0, mHaplotypeScorer.Score(ref_hap)));
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

// NOLINTNEXTLINE(readability-function-size)  // TODO(lancet): refactor to reduce function size
void SequenceComplexityScorer::ScoreDeltas(SequenceComplexity& cplx, std::string_view ref_hap,
                                           usize const ref_pos, usize const ref_len,
                                           std::string_view alt_hap, usize const alt_pos,
                                           usize const alt_len) const {
  // HRun delta at ±5bp
  auto const ref_hrun_window = ExtractFlank(ref_hap, ref_pos, ref_len, DELTA_HRUN_FLANK);
  auto const alt_hrun_window = ExtractFlank(alt_hap, alt_pos, alt_len, DELTA_HRUN_FLANK);
  cplx.mDeltaHRun = MaxHomopolymerRun(alt_hrun_window) - MaxHomopolymerRun(ref_hrun_window);

  // Entropy delta at ±10bp
  auto const ref_ent_window = ExtractFlank(ref_hap, ref_pos, ref_len, DELTA_ENTROPY_FLANK);
  auto const alt_ent_window = ExtractFlank(alt_hap, alt_pos, alt_len, DELTA_ENTROPY_FLANK);
  cplx.mDeltaEntropy = LocalShannonEntropy(alt_ent_window) - LocalShannonEntropy(ref_ent_window);

  // LongdustQ delta at ±50bp (log-space)
  auto const alt_lq_window = ExtractFlank(alt_hap, alt_pos, alt_len, LQ_FLANK);
  auto const alt_lq = std::log1p(std::max(0.0, mFlankScorer.Score(alt_lq_window)));
  cplx.mDeltaFlankLQ = alt_lq - cplx.mContextFlankLQ;
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

void SequenceComplexityScorer::ScoreTrMotif(SequenceComplexity& cplx, std::string_view alt_hap,
                                            usize const alt_pos, usize const alt_len) {
  auto const window = ExtractFlank(alt_hap, alt_pos, alt_len, TR_MOTIF_FLANK);
  auto const start = std::max(static_cast<i64>(0), static_cast<i64>(alt_pos) - TR_MOTIF_FLANK);
  auto const var_pos_in_window = static_cast<i32>(static_cast<i64>(alt_pos) - start);

  VariantTRFeatures trep;
  AccumulateTRFeatures(trep, window, var_pos_in_window, static_cast<i32>(alt_len),
                       static_cast<i32>(window.size()));

  // Sentinel-safe transforms: dist < 0 means no TR found → all zeros
  if (trep.mDistToNearestTr < 0) {
    cplx.mTrAffinity = 0.0F;
    cplx.mTrPurity = 0.0F;
    cplx.mTrPeriod = 0;
  } else {
    cplx.mTrAffinity = 1.0F / (1.0F + static_cast<f32>(trep.mDistToNearestTr));
    cplx.mTrPurity = trep.mNearestTrPurity;
    cplx.mTrPeriod = trep.mNearestTrPeriod;
  }
  cplx.mIsStutterIndel = trep.mIsStutterIndel;
}

// ============================================================================
// SequenceComplexity::FormatVcfValue — 11 comma-separated values
//
// Order: Context(4), Delta(3), TrMotif(4)
// ============================================================================

auto SequenceComplexity::FormatVcfValue() const -> std::string {
  return absl::StrCat(mContextHRun, ",", FormatComplexityScore(mContextEntropy), ",",
                      FormatComplexityScore(mContextFlankLQ), ",",
                      FormatComplexityScore(mContextHaplotypeLQ), ",", mDeltaHRun, ",",
                      FormatComplexityScore(mDeltaEntropy), ",",
                      FormatComplexityScore(mDeltaFlankLQ), ",", FormatComplexityScore(mTrAffinity),
                      ",", FormatComplexityScore(mTrPurity), ",", mTrPeriod, ",", mIsStutterIndel);
}

// ============================================================================
// SequenceComplexity::MergeMax — element-wise max for multi-haplotype
// and multi-allelic merging.
//
// Context features: max gives highest brittleness.
// Delta features: max gives most extreme perturbation.
// TR motif features: max gives closest/purest repeat.
// ============================================================================

void SequenceComplexity::MergeMax(SequenceComplexity const& other) {
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
