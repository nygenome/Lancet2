#include "lancet/caller/genotyper.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <utility>

extern "C" {
#include "mmpriv.h"
}

#include "absl/hash/hash.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_set.h"
#include "lancet/hts/phred_quality.h"


namespace {
  
  // ── Scoring constants for Illumina read-to-contig realignment ────────────────
  static constexpr int SCORING_MATCH = 1;
  static constexpr int SCORING_MISMATCH = 4;
  static constexpr int SCORING_GAP_OPEN = 12;
  static constexpr int SCORING_GAP_EXTEND = 3;
  static constexpr i8 ALPHABET_SIZE = 5;
  
  // 5×5 scoring matrix for ComputeLocalScore: A=0, C=1, G=2, T=3, N=4
  static constexpr auto MakeScoringMatrix() -> std::array<i8, 25> {
    std::array<i8, 25> mat{};
    for (i8 i = 0; i < ALPHABET_SIZE; ++i) {
      for (i8 j = 0; j < ALPHABET_SIZE; ++j) {
        if (i == ALPHABET_SIZE - 1 || j == ALPHABET_SIZE - 1) {
          mat[i * ALPHABET_SIZE + j] = 0;
        } else {
          mat[i * ALPHABET_SIZE + j] = (i == j)
              ? static_cast<i8>(SCORING_MATCH)
              : static_cast<i8>(-SCORING_MISMATCH);
        }
      }
    }
    return mat;
  }
  
  static constexpr std::array<i8, 25> SCORING_MATRIX = MakeScoringMatrix();
  
  // ASCII → numeric base encoding: A/a→0, C/c→1, G/g→2, T/t→3, else→4 (N)
  static constexpr auto MakeEncodingTable() -> std::array<u8, 256> {
    std::array<u8, 256> tbl{};
    for (auto& val : tbl) {
      val = 4;
    }
    tbl['A'] = 0; tbl['a'] = 0;
    tbl['C'] = 1; tbl['c'] = 1;
    tbl['G'] = 2; tbl['g'] = 2;
    tbl['T'] = 3; tbl['t'] = 3;
    return tbl;
  }
  
  static constexpr std::array<u8, 256> ENCODE_TABLE = MakeEncodingTable();

  // Free minimap2 alignment results (mm_reg1_t array)
  inline void FreeMm2Alignment(mm_reg1_t* regs, const int num_regs) {
    if (regs == nullptr) return;
    for (int idx = 0; idx < num_regs; ++idx) {
      std::free(regs[idx].p);
    }
    std::free(regs);
  }
  
}  // namespace

// ============================================================================
// ComputeLocalScore: evaluate alignment quality in a variant's region.
//
// Given a read→haplotype CIGAR alignment, this function extracts specific metrics
// for the sub-region of the haplotype that contains the variant natively:
//
//   1. pbq_score:  PBQ-weighted DP score. Each position's substitution matrix
//                  contribution is scaled by (1 - ε) where ε = 10^(-PBQ/10).
//                  This down-weights low-confidence bases and up-weights
//                  high-confidence ones, analogous to GATK's PairHMM which
//                  bakes PBQ directly into the per-read log-likelihood.
//
//   2. identity:   fraction of aligned bases that are exact matches.
//
//   4. base_qual:  Minimum Phred base quality (weakest-link) or flanking boundary.
//
// CRITICAL: tpos coordinates in the CIGAR are relative to the alignment start
// (ref_start from mm_map), NOT position 0 of the haplotype. The caller must
// adjust var_start by subtracting ref_start before calling this function.
//
//   Haplotype:   |----[var_start..........var_end)------|
//                      ^ref_start
//   CIGAR tpos:  0  1  2 ...
//                      ^var_start_in_aln = var_start - ref_start
// ============================================================================
namespace {

/// Convert a Phred quality score to confidence weight: 1 - 10^(-Q/10).
/// Q=0 → 0.0, Q=10 → 0.9, Q=20 → 0.99, Q=30 → 0.999, Q=40 → 0.9999
inline auto PhredToConfidence(const u8 qual) -> f64 {
  return 1.0 - lancet::hts::PhredToErrorProb(qual);
}

struct LocalScoreResult {
  f64 pbq_score = 0.0;
  f64 raw_score = 0.0; // Captures identical unweighted matrix paths identically
  f64 identity = 0.0;
  u8 base_qual = 0;
};

// ============================================================================
// RegionAccumulator: Local variant scoring abstraction state.
//
// Encapsulates scoring so the CIGAR walk loop stays clean. Evaluates alignment 
// quality specifically within a variant's physical array boundaries via absolute 
// haplotypic mapping rather than relative offsets.
// ============================================================================
struct RegionAccumulator {
  const absl::Span<const u8> mQuery;           // 16B
  const absl::Span<const u8> mTarget;          // 16B
  const absl::Span<const u8> mBaseQuals;       // 16B
  const std::array<i8, 25>& mScoringMatrix;    // 8B
  f64 mPbqScore = 0.0;                         // 8B
  f64 mRawScore = 0.0;                         // 8B
  usize mMatches = 0;                          // 8B
  usize mAligned = 0;                          // 8B
  const i32 mAlnRefStart;                      // 4B
  const i32 mVarStartHap;                      // 4B
  const i32 mVarEndHap;                        // 4B
  u8 mMinBq = 255;                             // 1B


  // Checks if the current alignment position overlaps the variant being scored.
  // 
  //   Haplotype Array :  [ A  B  C  D  E  F  G  H  I ]
  //   Variant Region  :           [var_start ... var_end)
  //   Alignment       :     [aln_start ... tpos_rel ... ]
  // 
  // We translate `tpos_rel` (which starts at 0 for the alignment) into an
  // absolute position (`abs_pos`) on the haplotype to see if it falls inside
  // the [mVarStartHap, mVarEndHap) window.
  [[nodiscard]] inline auto InRegion(i32 tpos_rel) const -> bool {
    const i32 abs_pos = mAlnRefStart + tpos_rel;
    return abs_pos >= mVarStartHap && abs_pos < mVarEndHap;
  }

  // Scores a single query-target base pair and adds it to the running totals.
  //   1. Raw Score: The penalty from the substitution matrix (e.g. mismatch = -4).
  //   2. PBQ Score: The raw penalty scaled by the base quality confidence.
  //                 A low quality mismatch is penalized less than a high quality one.
  //   3. Identity : Tracks pure exact matches to compute the exact-match fraction.
  void ScoreAlignedPair(i32 tpos_rel, usize qpos) {
    ++mAligned;
    if (qpos >= mQuery.size() || static_cast<usize>(tpos_rel) >= mTarget.size()) return;

    const auto raw = mScoringMatrix[mTarget[tpos_rel] * 5 + mQuery[qpos]];
    mRawScore += static_cast<f64>(raw);
    const f64 weight = (qpos < mBaseQuals.size()) ? PhredToConfidence(mBaseQuals[qpos]) : 1.0;
    mPbqScore += static_cast<f64>(raw) * weight;
    mMatches += static_cast<usize>(mQuery[qpos] == mTarget[tpos_rel]);
  }

  // Track the minimum base quality at an isolated read position natively.
  inline void TrackBaseQual(usize qpos) {
    if (qpos < mBaseQuals.size()) mMinBq = std::min(mMinBq, mBaseQuals[qpos]);
  }

  // Deletions do not have their own quality scores since the bases are missing.
  // Instead, we estimate the deletion's confidence by checking the quality of
  // the adjacent bases immediately surrounding the deletion.
  inline void TrackDeletionBounds(usize qpos) {
    if (qpos > 0 && qpos - 1 < mBaseQuals.size()) mMinBq = std::min(mMinBq, mBaseQuals[qpos - 1]);
    if (qpos < mBaseQuals.size()) mMinBq = std::min(mMinBq, mBaseQuals[qpos]);
  }

  // Finalize stats. min_bq defaults to 0 if we lacked base quality evidence entirely.
  [[nodiscard]] auto ToResult() const -> LocalScoreResult {
    return {
        .pbq_score = mPbqScore,
        .raw_score = mRawScore,
        .identity = mAligned > 0 ? static_cast<f64>(mMatches) / static_cast<f64>(mAligned) : 0.0,
        .base_qual = (mMinBq == 255) ? u8{0} : mMinBq,
    };
  }
};

auto ComputeLocalScore(const std::vector<lancet::hts::CigarUnit>& qry_aln_cigar,
                       absl::Span<const u8> qry_seq, absl::Span<const u8> hap_seq,
                       absl::Span<const u8> qry_quals,
                       i32 aln_start_on_hap, i32 var_start_on_hap, i32 var_len_on_hap,
                       const std::array<i8, 25>& score_matrix) -> LocalScoreResult {
  if (qry_aln_cigar.empty() || var_len_on_hap == 0) return {};

  RegionAccumulator acc{
      .mQuery = qry_seq,
      .mTarget = hap_seq,
      .mBaseQuals = qry_quals,
      .mScoringMatrix = score_matrix,
      .mAlnRefStart = aln_start_on_hap,
      .mVarStartHap = var_start_on_hap,
      .mVarEndHap = var_start_on_hap + var_len_on_hap,
  };
  i32 tpos = 0;
  usize qpos = 0;

  for (const auto& unit : qry_aln_cigar) {
    if (aln_start_on_hap + tpos >= acc.mVarEndHap && unit.ConsumesReference()) break;

    const auto op = unit.Operation();
    const u32 len = unit.Length();

    switch (op) {
      // ── Substitution Mappings ──
      // Consumes both query and target bases. Iteratively maps bases through the 
      // strict raw substitution matrix DP to calculate precise alignment identity 
      // and mismatch penalties explicitly restricted inside the variant boundary.
      case lancet::hts::CigarOp::ALIGNMENT_MATCH:
      case lancet::hts::CigarOp::SEQUENCE_MATCH:
      case lancet::hts::CigarOp::SEQUENCE_MISMATCH: {
        for (u32 i = 0; i < len; ++i, ++tpos, ++qpos) {
          if (!acc.InRegion(tpos)) continue;
          acc.ScoreAlignedPair(tpos, qpos);
          acc.TrackBaseQual(qpos);
        }
        break;
      }

      // ── Insertion Geometry ──
      // Consumes solely query bases. Penalizes the alignment quality natively via 
      // gap extension weights because these isolated query nucleotides lack reference backing.
      case lancet::hts::CigarOp::INSERTION: {
        // Inserted bases don't advance tpos, but if we're inside the variant
        // region they count as aligned content and contribute PBQ.
        const bool in_region = acc.InRegion(tpos);
        for (u32 i = 0; i < len; ++i, ++qpos) {
          if (!in_region) continue;
          ++acc.mAligned;
          acc.TrackBaseQual(qpos);
          acc.mRawScore += static_cast<f64>(SCORING_GAP_EXTEND);
          acc.mPbqScore += static_cast<f64>(SCORING_GAP_EXTEND);
        }
        break;
      }

      // ── Deletion Geometry ──
      // Consumes solely target bases. Symmetrically penalizes via gap extensions.
      // Since true deletions physically lack query nucleotides mathematically, the quality 
      // score relies explicitly on borrowing confidence from the structurally flanking neighbors.
      case lancet::hts::CigarOp::DELETION: {
        for (u32 i = 0; i < len; ++i, ++tpos) {
          if (acc.InRegion(tpos)) {
              ++acc.mAligned;
              acc.mRawScore += static_cast<f64>(SCORING_GAP_EXTEND);
              acc.mPbqScore += static_cast<f64>(SCORING_GAP_EXTEND);
          }
        }
        acc.TrackDeletionBounds(qpos);
        break;
      }

      // ── Unmapped Sequences ──
      // Adjusts spatial limits flawlessly without DP mutation scoring interactions.
      // Soft-clip toxicity defaults to being handled explicitly inside ComputeSoftClipPenalty.
      case lancet::hts::CigarOp::SOFT_CLIP:      { qpos += len; break; }
      case lancet::hts::CigarOp::REFERENCE_SKIP: { tpos += len; break; }
      default: break;
    }
  }

  return acc.ToResult();
}

// Calculates the penalty for unaligned sequence at the ends of a read.
// 
//   Read:  [Soft Clip] ====== Aligned Sequence ====== [Soft Clip]
//          5' End                                     3' End
// 
// We treat every soft-clipped base as a mismatch to penalize
// reads that only partially align to the haplotype.
inline auto ComputeSoftClipPenalty(const std::vector<lancet::hts::CigarUnit>& cigar) -> f64 {
  if (cigar.empty()) return 0.0;

  const auto& first = cigar.front();
  const auto& last = cigar.back();

  // Branchless boolean arithmetic: true translates to 1, false translates to 0. 
  // Prevents CPU pipeline stalls by avoiding conditional jump instructions natively.
  const i32 unaligned_5p = (first.Operation() == lancet::hts::CigarOp::SOFT_CLIP) * first.Length();
  const i32 unaligned_3p = (cigar.size() > 1) * (last.Operation() == lancet::hts::CigarOp::SOFT_CLIP) * last.Length();

  return static_cast<f64>(unaligned_5p + unaligned_3p) * SCORING_MISMATCH;
}

}  // namespace


namespace lancet::caller {

// ============================================================================
// Constructor: initialize minimap2 with custom Illumina scoring parameters.
//
// We do NOT use the 'sr' preset because its gap penalties are tuned for
// whole-genome alignment where gaps represent biological variation. In local
// assembly, biological variation is already in the haplotype sequences —
// gaps here are machine errors and must be heavily penalized.
//
// See the SCORING_* constants defined natively in the anonymous namespace above for the full rationale.
// ============================================================================
Genotyper::Genotyper() {
  // 0 -> no info, 1 -> error, 2 -> warning, 3 -> debug
  mm_verbose = 1;

  // Start from the default parameter set, then override scoring
  mm_set_opt(nullptr, mIndexingOpts.get(), mMappingOpts.get());

  auto* mopts = mMappingOpts.get();
  mopts->a = SCORING_MATCH;
  mopts->b = SCORING_MISMATCH;
  mopts->q = SCORING_GAP_OPEN;
  mopts->e = SCORING_GAP_EXTEND;
  // Use single-affine gap model: set the second model to the same parameters.
  // This disables minimap2's dual-affine (convex) model.
  mopts->q2 = SCORING_GAP_OPEN;
  mopts->e2 = SCORING_GAP_EXTEND;

  // 1. Z-Drop (zdrop / zdrop_inv)
  // Minimap2 terminates DP alignment if the local Smith-Waterman score drops below
  // (max_score - zdrop). If a tumor contains a massive 300bp somatic deletion,
  // the affine gap penalty (e.g. O + 300*E) structurally exceeds standard Z-drop thresholds
  // (default 400). This causes minimap2 to silently truncate the alignment and report
  // separate supplementary reads rather than a contiguous CIGAR spanning the deletion!
  // By overriding Z-drop to 100000, we mathematically disable alignment truncation,
  // forcing the band alignment to traverse any arbitrarily large intra-window structural variations.
  mopts->zdrop = 100000;
  mopts->zdrop_inv = 100000;

  // 2. Bandwidth (bw)
  // Minimap2 only computes DP matrix bounds within `bw` distance from the main diagonal
  // (O(N * bw) complexity). If an insertion is larger than `bw`, the alignment mathematically
  // strikes the banding boundary and terminates or forces chaotic mismatches. Because Lancet2's
  // active regions encompass massive micro-assemblies (where germline insertions can exceed 2kb), 
  // we must drastically inflate `bw` from its default (often 500) to 10000. This guarantees 
  // massive germline and somatic insertions never exceed the matrix diagonal thresholds.
  mopts->bw = 10000;

  // 3. Seeding Sensitivity (k, w)
  // Minimap2 uses (k, w)-minimizers to detect identical seed anchors before executing SW.
  // In highly mutated sequences, STRs, or clustered hypermutation sites, the default 
  // initialization (k=15, w=10) drops anchors because a continuous 15bp exact match 
  // rarely exists. By overriding this down to k=11 and w=5, we dramatically increase 
  // algorithmic sensitivity, natively instantiating index anchors directly through 
  // deeply fragmented micro-windows.
  mIndexingOpts->k = 11;
  mIndexingOpts->w = 5;

  mopts->flag |= MM_F_CIGAR;  // Generate CIGAR (needed for local scoring)
  mopts->best_n = 1;          // Only keep the best hit per haplotype
}

// ============================================================================
// Genotype: main entry point — orchestrates alignment and evidence collection.
//
//  ┌──────────────┐    ┌──────────────┐    ┌─────────────┐
//  │  Haplotypes  │    │    Reads     │    │  VariantSet │
//  │ (REF + ALTs) │    │  (all samps) │    │ (raw vars)  │
//  └───────┬──────┘    └──────┬───────┘    └─────┬───────┘
//          │                  │                  │
//          ▼                  │                  │
//   ResetData()               │                  │
//   (build mm2 indices)       │                  │
//          │                  │                  │
//          │    ┌─────────────┘                  │
//          │    │                                │
//          ▼    ▼                                │
//   AssignReadToAlleles() ◄──────────────────────┘
//   (mm_map per hap,
//    local scoring per var)
//          │
//          ▼
//   AddToTable()
//   (ReadEvidence → VariantSupport)
//          │
//          ▼
//   ┌──────────────┐
//   │    Result    │
//   │  per-variant │
//   │  per-sample  │
//   │  support     │
//   └──────────────┘
// ============================================================================
auto Genotyper::Genotype(Haplotypes hap_seqs, Reads qry_reads, const VariantSet& variant_set) -> Result {
  ResetData(hap_seqs);

  Result out_vars_table;
  for (const auto& qry_read : qry_reads) {
    auto allele_assignments = AssignReadToAlleles(qry_read, variant_set);
    AddToTable(out_vars_table, qry_read, allele_assignments);
  }

  return out_vars_table;
}


// ============================================================================
// ResetData: build minimap2 indices for all haplotype sequences.
//
// Each haplotype gets its own index so we can align reads independently
// to REF and each ALT haplotype and compare scores.
// ============================================================================
void Genotyper::ResetData(Haplotypes hap_seqs) {
  mIndices.clear();
  mIndices.reserve(hap_seqs.size());

  const auto* iopts = mIndexingOpts.get();
  for (const auto& hap_seq : hap_seqs) {
    const char* raw_seq = hap_seq.c_str();
    auto* idx_result = mm_idx_str(iopts->w, iopts->k, 0, iopts->bucket_bits, 1, &raw_seq, nullptr);
    mIndices.emplace_back(Minimap2Index(idx_result));
  }

  // Pre-encode haplotype sequences for local scoring.
  // mm_idx stores sequences internally but doesn't expose them via a clean API,
  // so we maintain our own numeric-encoded copies for ComputeLocalScore.
  mEncodedHaplotypes.clear();
  mEncodedHaplotypes.reserve(hap_seqs.size());
  for (const auto& hap_seq : hap_seqs) {
    mEncodedHaplotypes.push_back(EncodeSequence(hap_seq));
  }

  auto* mopts = mMappingOpts.get();
  for (const auto& mm2_idx : mIndices) {
    mm_mapopt_update(mopts, mm2_idx.get());
  }
}

auto Genotyper::AssignReadToAlleles(const cbdg::Read& qry_read,
                                    const VariantSet& variant_set) -> PerVariantAssignment {
  auto all_alns = AlignToAllHaplotypes(qry_read);
  if (all_alns.empty()) return {};

  const auto qry_quals = qry_read.QualView();
  const auto qry_seq_encoded = EncodeSequence(qry_read.SeqView());
  const usize qry_read_length = qry_read.Length();
  PerVariantAssignment allele_assignments;

  // O(N) PERFORMANCE WIN: Extracted from the variant iterator loop.
  // Calculated natively exactly once per read!
  const u32 baseline_ref_nm = ComputeRefEditDistance(all_alns, absl::MakeConstSpan(qry_seq_encoded), qry_read_length);

  for (const auto& variant : variant_set) {
    ReadAlleleAssignment best{};
    best.ref_nm = baseline_ref_nm;
    best.allele = REF_ALLELE_IDX;
    best.global_score = std::numeric_limits<i32>::lowest();
    
    f64 best_combined = std::numeric_limits<f64>::lowest();
    bool found_any = false;

    // Scan the read natively against all identically aligned physical haplotypes
    for (const auto& aln : all_alns) {
      found_any |= EvaluateAlignment(aln, variant, absl::MakeConstSpan(qry_seq_encoded), qry_quals, 
                                     qry_read_length, best, best_combined);
    }

    if (found_any) {
      allele_assignments.emplace(&variant, best);
    }
  }

  return allele_assignments;
}

// ============================================================================
// AlignToAllHaplotypes: Exhaustively map a read to all haplotype sequences natively.
//
// Uses minimap2's full pipeline (seed → chain → extend) rather than raw ksw2 DP 
// because minimizer hashing efficiently resolves topological boundaries natively before 
// falling back onto dynamic programming execution.
//
// CRITICAL: We explicitly do NOT implement an early-exit short-circuit (even if a 
// read perfectly matches a haplotype). We must deliberately exhaust the entire 
// Haplotype Space to mathematically guarantee we assemble valid cross-haplotype noise 
// constraints, compute relative Edit Distances mapped strictly against the true Reference 
// baseline natively, and properly enforce global best-match boundaries unconditionally.
// ============================================================================
auto Genotyper::AlignToAllHaplotypes(const cbdg::Read& qry_read) -> std::vector<Mm2AlnResult> {
  std::vector<Mm2AlnResult> results;
  results.reserve(mIndices.size());

  int nregs = 0;
  auto* tbuffer = mThreadBuffer.get();
  const auto* map_opts = mMappingOpts.get();
  const auto read_len = static_cast<int>(qry_read.Length());

  for (usize idx = 0; idx < mIndices.size(); ++idx) {
    const auto* hap_mm_idx = mIndices[idx].get();
    auto* regs = mm_map(hap_mm_idx, read_len, qry_read.SeqPtr(), &nregs, tbuffer, map_opts, qry_read.QnamePtr());

    if (regs == nullptr || nregs <= 0) {
      FreeMm2Alignment(regs, nregs);
      continue;
    }

    // Take the top hit only (best_n = 1)
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const mm_reg1_t* top_hit = &regs[0];

    Mm2AlnResult result;
    result.score = top_hit->score;
    result.ref_start = top_hit->rs;  // critical: where alignment starts on haplotype
    result.ref_end = top_hit->re;
    result.identity = mm_event_identity(top_hit);
    result.hap_idx = idx;

    // Extract CIGAR from mm_extra_t
    if (top_hit->p != nullptr && top_hit->p->n_cigar > 0) {
      result.cigar.reserve(top_hit->p->n_cigar);
      for (u32 k = 0; k < top_hit->p->n_cigar; ++k) {
        result.cigar.emplace_back(top_hit->p->cigar[k]);
      }
    }

    results.push_back(std::move(result));
    FreeMm2Alignment(regs, nregs);
  }

  return results;
}

// ============================================================================
// EncodeSequence: ASCII DNA → numeric (0-4) for local scoring.
// Single pass using the constexpr lookup table. O(n), no branches.
// ============================================================================
auto Genotyper::EncodeSequence(const std::string_view raw_seq) -> std::vector<u8> {
  std::vector<u8> encoded(raw_seq.size());
  for (usize i = 0; i < raw_seq.size(); ++i) {
    encoded[i] = ENCODE_TABLE[static_cast<u8>(raw_seq[i])];
  }
  return encoded;
}

auto Genotyper::ComputeRefEditDistance(const std::vector<Mm2AlnResult>& alns,
                                       absl::Span<const u8> qry_seq_encoded, 
                                       usize qry_read_length) const -> u32 {
  // Calculate the edit distance (NM) against the Reference Haplotype.
  // We explicitly compute this against the REF sequence, even if the read
  // was assigned to an ALT allele. This provides a baseline noise metric 
  // used later to filter out false positives.
  for (const auto& aln : alns) {
    if (aln.hap_idx != REF_HAP_IDX || aln.ref_start >= aln.ref_end) continue;
    const auto ref_aln_len = static_cast<usize>(aln.ref_end - aln.ref_start);
    const auto ref_target = absl::MakeConstSpan(mEncodedHaplotypes[REF_HAP_IDX])
                                .subspan(static_cast<usize>(aln.ref_start), ref_aln_len);
    u32 nm = hts::ComputeEditDistance(aln.cigar, qry_seq_encoded, ref_target);
    
    for (const auto& u : aln.cigar) {
      if (u.Operation() == hts::CigarOp::SOFT_CLIP) nm += u.Length();
    }
    return nm;
  }
  // Initialize the edit distance (ref_nm) to the maximum possible value 
  // (the full read length). If we find a valid alignment later, we'll 
  // replace this with the real edit distance.
  return static_cast<u32>(qry_read_length);
}

auto Genotyper::EvaluateAlignment(const Mm2AlnResult& aln, 
                                  const RawVariant& variant,
                                  absl::Span<const u8> qry_seq_encoded, 
                                  absl::Span<const u8> qry_base_quals, 
                                  usize qry_read_length,
                                  ReadAlleleAssignment& best, 
                                  f64& best_combined) const -> bool {
  i32 hap_var_start = 0, hap_var_len = 0;
  AlleleIndex mapped_allele_idx = REF_ALLELE_IDX;

  if (!ExtractHapBounds(variant, aln.hap_idx, hap_var_start, hap_var_len, mapped_allele_idx)) return false;

  const i32 hap_var_end = hap_var_start + hap_var_len;
  // Skip this alignment if it doesn't physically overlap the variant region.
  //
  //   Alignment Span :            [ref_start ...... ref_end)
  //   Variant A      :  [start..end)                           (Skips)
  //   Variant B      :                 [start..end)            (Processes)
  //   Variant C      :                                  [start..end) (Skips)
  if (hap_var_end <= aln.ref_start || hap_var_start >= aln.ref_end) return false;

  const auto aln_len = static_cast<usize>(aln.ref_end - aln.ref_start);
  const auto target = absl::MakeConstSpan(mEncodedHaplotypes[aln.hap_idx])
                          .subspan(static_cast<usize>(aln.ref_start), aln_len);

  const auto local = ComputeLocalScore(
      aln.cigar, qry_seq_encoded, target,
      qry_base_quals, aln.ref_start, hap_var_start, hap_var_len, SCORING_MATRIX);

  // Subtract the soft-clip penalty from the global score. 
  // This prevents supplementary alignments (where half the read is soft-clipped)
  // from incorrectly outscoring full contiguous alignments.
  const f64 sc_penalty = ComputeSoftClipPenalty(aln.cigar);
  const f64 global_adjusted = static_cast<f64>(aln.score) - sc_penalty;
  const f64 combined = global_adjusted - local.raw_score + (local.pbq_score * local.identity);

  if (best_combined == std::numeric_limits<f64>::lowest() || combined > best_combined) {
    best.allele = mapped_allele_idx;
    best.global_score = static_cast<i32>(global_adjusted - local.raw_score);
    best.local_score = local.pbq_score;
    best.local_identity = local.identity;
    best.base_qual_at_var = local.base_qual;
    
    // ── Folded read position ──────────────────────────────────
    usize var_start_in_aln = 0;
    if (hap_var_start > aln.ref_start) var_start_in_aln = static_cast<usize>(hap_var_start - aln.ref_start);
    
    const auto qpos_at_var = hts::CigarRefPosToQueryPos(aln.cigar, var_start_in_aln);
    const auto rel_pos = qry_read_length > 0 ? static_cast<f64>(qpos_at_var) / static_cast<f64>(qry_read_length) : 0.5;
    best.folded_read_pos = std::min(rel_pos, 1.0 - rel_pos);

    best_combined = combined;
    return true;
  }

  return false;
}

// ============================================================================
// AssignReadToAlleles
//
// Determines which allele (REF or ALT) a read best supports for each variant.
// Evaluates the read against all possible haplotypes and scores the fit using:
//   Score = (Global_Score - Local_Raw) + (Local_PBQ * Local_Identity)
// ============================================================================
// ────────────────────────────────────────────────────────────────────────────
// AssignReadToAlleles Internal Core Helpers
// ────────────────────────────────────────────────────────────────────────────
auto Genotyper::ExtractHapBounds(const RawVariant& variant, usize aln_hap_idx,
                                 i32& out_hap_variant_start, i32& out_hap_variant_length, AlleleIndex& out_allele) const -> bool {
  if (aln_hap_idx == REF_HAP_IDX) {
    out_hap_variant_start = static_cast<i32>(variant.mLocalRefStart0Idx);
    out_hap_variant_length = static_cast<i32>(variant.mRefAllele.size());
    out_allele = REF_ALLELE_IDX;
    return true;
  }
  for (usize alt_pos = 0; alt_pos < variant.mAlts.size(); ++alt_pos) {
    const auto& alt_allele = variant.mAlts[alt_pos];
    auto it = alt_allele.mLocalHapStart0Idxs.find(aln_hap_idx);
    if (it != alt_allele.mLocalHapStart0Idxs.end()) {
      out_hap_variant_start = static_cast<i32>(it->second);
      // We use the exact length of the ALT allele sequence because
      // multiallelics can have different lengths than the REF allele.
      out_hap_variant_length = static_cast<i32>(alt_allele.mSequence.size());
      out_allele = AlleleIndex(alt_pos + 1);
      return true;
    }
  }
  return false;
}

// ============================================================================
// AddToTable: record a read's allele assignments into per-variant support.
//
// Constructs ReadEvidence with all available read-level metrics:
//   - base quality at the variant position
//   - original mapping quality
//   - normalized alignment score
//   - strand direction
// All metrics flow through to VCF FORMAT fields via VariantSupport aggregation.
// ============================================================================
void Genotyper::AddToTable(Result& out_vars_table, const cbdg::Read& qry_read,
                           const PerVariantAssignment& allele_assignments) {
  const auto sample_name = qry_read.SampleName();
  const auto rname_hash = absl::HashOf(qry_read.QnameView());
  const auto strand = qry_read.BitwiseFlag().IsRevStrand() ? Strand::REV : Strand::FWD;

  for (const auto& [var_ptr, assignment] : allele_assignments) {
    // Look up (or create) the per-sample evidence aggregator for this variant.
    // Result is keyed: variant → sample_name → VariantSupport.
    // Default insertion occurs seamlessly at both tier levels:
    //   rslt[var_ptr]               → inserts an empty SupportArray interface
    //   .FindOrCreate(sample_name)  → intercepts key searches, appending and allocating 
    //                                 a unique_ptr<VariantSupport> on-the-fly dynamically.
    auto& support = out_vars_table[var_ptr].FindOrCreate(sample_name);

    VariantSupport::ReadEvidence evidence;
    evidence.rname_hash = static_cast<u32>(rname_hash);
    evidence.allele = assignment.allele;
    evidence.strand = strand;
    evidence.base_qual = assignment.base_qual_at_var;
    evidence.map_qual = qry_read.MapQual();
    evidence.aln_score = static_cast<f64>(assignment.CombinedScore());
    evidence.folded_read_pos = assignment.folded_read_pos;
    evidence.is_soft_clipped = qry_read.IsSoftClipped();
    evidence.is_proper_pair = qry_read.IsProperPair();
    evidence.insert_size = qry_read.InsertSize();
    evidence.ref_nm = assignment.ref_nm;

    support.AddEvidence(evidence);
  }
}

}  // namespace lancet::caller
