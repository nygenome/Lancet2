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

namespace lancet::caller {

// ============================================================================
// Constructor: initialize minimap2 with custom Illumina scoring parameters.
//
// We do NOT use the 'sr' preset because its gap penalties are tuned for
// whole-genome alignment where gaps represent biological variation. In local
// assembly, biological variation is already in the haplotype sequences —
// gaps here are machine errors and must be heavily penalized.
//
// See the SCORING_* constants in genotyper.h for the full rationale.
// ============================================================================
Genotyper::Genotyper() {
  // 0 -> no info, 1 -> error, 2 -> warning, 3 -> debug
  mm_verbose = 1;

  // Start from the default parameter set, then override scoring
  mm_set_opt(nullptr, mIndexingOpts.get(), mMappingOpts.get());

  auto* mopts = mMappingOpts.get();
  mopts->a = genotyper_detail::SCORING_MATCH;
  mopts->b = genotyper_detail::SCORING_MISMATCH;
  mopts->q = genotyper_detail::SCORING_GAP_OPEN;
  mopts->e = genotyper_detail::SCORING_GAP_EXTEND;
  // Use single-affine gap model: set the second model to the same parameters.
  // This disables minimap2's dual-affine (convex) model.
  mopts->q2 = genotyper_detail::SCORING_GAP_OPEN;
  mopts->e2 = genotyper_detail::SCORING_GAP_EXTEND;

  mopts->flag |= MM_F_CIGAR;  // Generate CIGAR (needed for local scoring)
  mopts->best_n = 1;          // Only keep the best hit per haplotype
}

// ============================================================================
// EncodeSequence: ASCII DNA → numeric (0-4) for local scoring.
// Single pass using the constexpr lookup table. O(n), no branches.
// ============================================================================
auto Genotyper::EncodeSequence(const std::string_view seq) -> std::vector<u8> {
  std::vector<u8> encoded(seq.size());
  for (usize i = 0; i < seq.size(); ++i) {
    encoded[i] = genotyper_detail::ENCODE_TABLE[static_cast<u8>(seq[i])];
  }
  return encoded;
}

// ============================================================================
// ResetData: build minimap2 indices for all haplotype sequences.
//
// Each haplotype gets its own index so we can align reads independently
// to REF and each ALT haplotype and compare scores.
// ============================================================================
void Genotyper::ResetData(Haplotypes sequences) {
  mIndices.clear();
  mIndices.reserve(sequences.size());

  const auto* iopts = mIndexingOpts.get();
  for (const auto& seq : sequences) {
    const char* raw_seq = seq.c_str();
    auto* idx_result = mm_idx_str(iopts->w, iopts->k, 0, iopts->bucket_bits, 1, &raw_seq, nullptr);
    mIndices.emplace_back(Minimap2Index(idx_result));
  }

  // Pre-encode haplotype sequences for local scoring.
  // mm_idx stores sequences internally but doesn't expose them via a clean API,
  // so we maintain our own numeric-encoded copies for ComputeLocalScore.
  mEncodedHaplotypes.clear();
  mEncodedHaplotypes.reserve(sequences.size());
  for (const auto& seq : sequences) {
    mEncodedHaplotypes.push_back(EncodeSequence(seq));
  }

  auto* mopts = mMappingOpts.get();
  for (const auto& mm2_idx : mIndices) {
    mm_mapopt_update(mopts, mm2_idx.get());
  }
}

namespace {

// Free minimap2 alignment results (mm_reg1_t array)
inline void FreeMm2Alignment(mm_reg1_t* regs, const int num_regs) {
  if (regs == nullptr) return;
  // NOLINTBEGIN(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc)
  for (int idx = 0; idx < num_regs; ++idx) {
    std::free(regs[idx].p);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  }
  std::free(regs);
  // NOLINTEND(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc)
}

}  // namespace

// ============================================================================
// AlignToAllHaplotypes: map a read to all haplotype indices via mm_map.
//
// Uses minimap2's full pipeline (seed → chain → extend) which is much faster
// than raw ksw2 DP for well-matching sequences because seeding resolves most
// of the alignment before DP is even invoked.
//
// Early-exit: if a read has a perfect match to any haplotype (identity ≈ 1.0,
// full query coverage), skip remaining haplotypes. This read unambiguously
// belongs to that haplotype and doesn't need comparison.
// ============================================================================
auto Genotyper::AlignToAllHaplotypes(const cbdg::Read& read) -> std::vector<Mm2AlnResult> {
  std::vector<Mm2AlnResult> results;
  results.reserve(mIndices.size());

  int nregs = 0;
  auto* tbuffer = mThreadBuffer.get();
  const auto* map_opts = mMappingOpts.get();
  const auto read_len = static_cast<int>(read.Length());

  for (usize idx = 0; idx < mIndices.size(); ++idx) {
    const auto* hap_mm_idx = mIndices[idx].get();
    auto* regs = mm_map(hap_mm_idx, read_len, read.SeqPtr(), &nregs, tbuffer, map_opts, read.QnamePtr());

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

    // Perfect match: every base of the read aligns to this haplotype with no mismatches.
    // qs==0 && qe==read_len means the full query is consumed (no soft clips).
    // identity==1.0 means zero mismatches/indels in the aligned region.
    result.is_full_match = (top_hit->qs == 0 && top_hit->qe == read_len && result.identity == 1.0);

    results.push_back(std::move(result));
    FreeMm2Alignment(regs, nregs);

    // Early exit: if perfect match to any haplotype, this read
    // unambiguously belongs here. Skip aligning to other haps.
    if (results.back().is_full_match) {
      break;
    }
  }

  return results;
}

// ============================================================================
// ComputeLocalScore: evaluate alignment quality in a variant's region.
//
// Given a read→haplotype CIGAR alignment, this function extracts three metrics
// for the sub-region of the haplotype that contains the variant:
//
//   1. score:      PBQ-weighted DP score. Each position's substitution matrix
//                  contribution is scaled by (1 - ε) where ε = 10^(-PBQ/10).
//                  This down-weights low-confidence bases and up-weights
//                  high-confidence ones, analogous to GATK's PairHMM which
//                  bakes PBQ directly into the per-read log-likelihood.
//
//   2. identity:   fraction of aligned bases that are exact matches.
//
//   3. base_qual:  minimum Phred base quality across all read positions
//                  that fall within the variant region (weakest-link).
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
  return 1.0 - std::pow(10.0, -static_cast<f64>(qual) / 10.0);
}

struct LocalScoreResult {
  f64 score = 0.0;      // PBQ-weighted DP score within variant region
  f64 identity = 0.0;   // matches / total_aligned in variant region
  u8 base_qual = 0;     // min PBQ across read positions in variant region
};

/// Accumulates alignment statistics within a variant region.
/// Encapsulates all scoring state so the CIGAR walk loop stays clean.
struct RegionAccumulator {
  const absl::Span<const u8> query;
  const absl::Span<const u8> target;
  const absl::Span<const u8> base_quals;
  const std::array<i8, 25>& scoring_matrix;
  const usize var_start;
  const usize var_end;

  f64 score = 0.0;
  usize matches = 0;
  usize aligned = 0;
  u8 min_bq = 255;

  [[nodiscard]] auto InRegion(usize tpos) const -> bool {
    return tpos >= var_start && tpos < var_end;
  }

  /// Score an aligned base pair at (tpos, qpos) within the variant region.
  void ScoreAlignedPair(usize tpos, usize qpos) {
    ++aligned;
    if (qpos >= query.size() || tpos >= target.size()) return;

    const auto raw = scoring_matrix[target[tpos] * 5 + query[qpos]];
    const f64 weight = (qpos < base_quals.size()) ? PhredToConfidence(base_quals[qpos]) : 1.0;
    score += static_cast<f64>(raw) * weight;
    matches += static_cast<usize>(query[qpos] == target[tpos]);
  }

  /// Track the minimum base quality at a query position.
  void TrackBaseQual(usize qpos) {
    if (qpos < base_quals.size()) {
      min_bq = std::min(min_bq, base_quals[qpos]);
    }
  }

  /// Finalize accumulated stats into a LocalScoreResult.
  [[nodiscard]] auto ToResult() const -> LocalScoreResult {
    return {
        .score = score,
        .identity = aligned > 0 ? static_cast<f64>(matches) / static_cast<f64>(aligned) : 0.0,
        // min_bq stays 255 if no query positions fell in the region (pure deletion).
        // Default to 0 = no base quality evidence available.
        // See ReadAlleleAssignment::base_qual_at_var in genotyper.h for rationale.
        .base_qual = (min_bq == 255) ? u8{0} : min_bq,
    };
  }
};

auto ComputeLocalScore(const std::vector<hts::CigarUnit>& cigar,
                       absl::Span<const u8> query, absl::Span<const u8> target,
                       absl::Span<const u8> base_quals,
                       const usize var_start, const usize var_len,
                       const std::array<i8, 25>& scoring_matrix) -> LocalScoreResult {
  if (cigar.empty() || var_len == 0) return {};

  RegionAccumulator acc{query, target, base_quals, scoring_matrix, var_start, var_start + var_len};
  usize tpos = 0;
  usize qpos = 0;

  for (const auto& unit : cigar) {
    if (tpos >= acc.var_end && unit.ConsumesReference()) break;

    const auto op = unit.Operation();
    const u32 len = unit.Length();

    switch (op) {
      case hts::CigarOp::ALIGNMENT_MATCH:
      case hts::CigarOp::SEQUENCE_MATCH:
      case hts::CigarOp::SEQUENCE_MISMATCH: {
        for (u32 i = 0; i < len; ++i, ++tpos, ++qpos) {
          if (!acc.InRegion(tpos)) continue;
          acc.ScoreAlignedPair(tpos, qpos);
          acc.TrackBaseQual(qpos);
        }
        break;
      }

      case hts::CigarOp::INSERTION: {
        // Inserted bases don't advance tpos, but if we're inside the variant
        // region they count as aligned content and contribute PBQ.
        const bool in_region = (tpos > var_start && tpos <= acc.var_end);
        for (u32 i = 0; i < len; ++i, ++qpos) {
          if (!in_region) continue;
          ++acc.aligned;
          acc.TrackBaseQual(qpos);
        }
        break;
      }

      case hts::CigarOp::DELETION: {
        // No query base → no PBQ contribution, just track aligned count.
        for (u32 i = 0; i < len; ++i, ++tpos) {
          if (acc.InRegion(tpos)) ++acc.aligned;
        }
        break;
      }

      case hts::CigarOp::SOFT_CLIP:      { qpos += len; break; }
      case hts::CigarOp::REFERENCE_SKIP:  { tpos += len; break; }
      default: break;
    }
  }

  return acc.ToResult();
}

}  // namespace

// ============================================================================
// AssignReadToAlleles: the Phase 2 abstraction boundary.
//
// For each variant in the VariantSet, determines which allele this read
// supports by comparing combined (global + local * identity) alignment scores
// across all haplotypes that carry that variant.
//
// CRITICAL FIX: mm_map's CIGAR is relative to ref_start (where the alignment
// begins on the haplotype), NOT position 0. So when computing local scores,
// the variant's position must be adjusted: var_start_in_aln = var_start - ref_start.
// If the variant falls outside the alignment span, it's skipped.
// ============================================================================
auto Genotyper::AssignReadToAlleles(const cbdg::Read& read,
                                    const VariantSet& vset) -> PerVariantAssignment {
  auto all_alns = AlignToAllHaplotypes(read);
  if (all_alns.empty()) {
    return {};
  }

  const auto base_quals = read.QualView();
  const auto read_seq = read.SeqView();
  const auto encoded_query = EncodeSequence(read_seq);

  PerVariantAssignment assignments;
  for (const auto& variant : vset) {
    ReadAlleleAssignment best{};
    best.allele = REF_ALLELE_IDX;
    best.global_score = 0;
    best.local_score = 0.0;
    best.local_identity = 0.0;
    best.base_qual_at_var = 0;
    bool found_any = false;

    for (const auto& aln : all_alns) {
      // Skip if this haplotype doesn't carry this variant
      const auto hap_pos_it = variant.mHapStart0Idxs.find(aln.hap_idx);
      if (hap_pos_it == variant.mHapStart0Idxs.end()) {
        continue;
      }

      const auto var_start_in_hap = static_cast<i32>(hap_pos_it->second);
      const auto var_len = std::max(variant.mRefAllele.size(), variant.mAltAllele.size());

      // CRITICAL: adjust variant position relative to alignment start.
      // mm_map's CIGAR starts at ref_start, not haplotype position 0.
      // If the variant falls outside the alignment span, skip.
      if (var_start_in_hap < aln.ref_start || var_start_in_hap >= aln.ref_end) {
        continue;
      }
      const auto var_start_in_aln = static_cast<usize>(var_start_in_hap - aln.ref_start);

      // Use pre-encoded haplotype sequence for local scoring.
      // The CIGAR covers [ref_start, ref_end) on the haplotype, so we pass
      // the corresponding subspan of the encoded haplotype as the target.
      const auto aln_len = static_cast<usize>(aln.ref_end - aln.ref_start);
      const auto target = absl::MakeConstSpan(mEncodedHaplotypes[aln.hap_idx])
                              .subspan(static_cast<usize>(aln.ref_start), aln_len);

      const auto local = ComputeLocalScore(
          aln.cigar, absl::MakeConstSpan(encoded_query), target,
          base_quals, var_start_in_aln, var_len, genotyper_detail::SCORING_MATRIX);

      // Use combined score (global + local * identity) for allele assignment.
      // See ReadAlleleAssignment comment block for rationale.
      const f64 combined = static_cast<f64>(aln.score) + (local.score * local.identity);
      if (!found_any || combined > best.CombinedScore()) {
        // Haplotype index ≠ allele index. Each RawVariant is bi-allelic:
        // REF haplotype (idx 0) → allele 0, any ALT haplotype → allele 1.
        // Multiple haplotypes may carry the same ALT allele.
        best.allele = (aln.hap_idx == REF_HAP_IDX) ? REF_ALLELE_IDX : AlleleIndex{1};
        best.global_score = aln.score;
        best.local_score = local.score;
        best.local_identity = local.identity;
        best.base_qual_at_var = local.base_qual;
        found_any = true;
      }
    }

    if (found_any) {
      assignments.emplace(&variant, best);
    }
  }

  return assignments;
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
void Genotyper::AddToTable(Result& rslt, const cbdg::Read& read,
                           const PerVariantAssignment& assignments) {
  const auto sample_name = read.SampleName();
  const auto rname_hash = absl::HashOf(read.QnameView());
  const auto strand = read.BitwiseFlag().IsRevStrand() ? Strand::REV : Strand::FWD;

  for (const auto& [var_ptr, assignment] : assignments) {
    // Look up (or create) the per-sample evidence aggregator for this variant.
    // Result is keyed: variant → sample_name → VariantSupport.
    // operator[] default-inserts at both levels if the key is absent:
    //   rslt[var_ptr]              → inserts empty PerSampleEvidence map
    //   rslt[var_ptr][sample_name] → inserts null unique_ptr<VariantSupport>
    // The null check below then constructs the VariantSupport on first access.
    auto& support = rslt[var_ptr][sample_name];
    if (!support) {
      support = std::make_unique<VariantSupport>();
    }

    VariantSupport::ReadEvidence evidence;
    evidence.rname_hash = static_cast<u32>(rname_hash);
    evidence.allele = assignment.allele;
    evidence.strand = strand;
    evidence.base_qual = assignment.base_qual_at_var;
    evidence.map_qual = read.MapQual();
    evidence.aln_score = static_cast<f64>(assignment.CombinedScore());

    support->AddEvidence(evidence);
  }
}

// ============================================================================
// Genotype: main entry point — orchestrates alignment and evidence collection.
//
//  ┌──────────────┐   ┌──────────────┐   ┌──────────────┐
//  │  Haplotypes   │   │    Reads     │   │  VariantSet  │
//  │ (REF + ALTs)  │   │  (all samps) │   │ (raw vars)   │
//  └───────┬───────┘   └──────┬───────┘   └──────┬───────┘
//          │                  │                   │
//          ▼                  │                   │
//   ResetData()               │                   │
//   (build mm2 indices)       │                   │
//          │                  │                   │
//          │    ┌─────────────┘                   │
//          │    │                                 │
//          ▼    ▼                                 │
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
//   │    Result     │
//   │  per-variant  │
//   │  per-sample   │
//   │  support      │
//   └──────────────┘
// ============================================================================
auto Genotyper::Genotype(Haplotypes haps, Reads reads, const VariantSet& vset) -> Result {
  ResetData(haps);

  Result result;
  for (const auto& read : reads) {
    auto assignments = AssignReadToAlleles(read, vset);
    AddToTable(result, read, assignments);
  }

  return result;
}

}  // namespace lancet::caller
