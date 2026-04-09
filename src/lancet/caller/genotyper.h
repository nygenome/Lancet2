#ifndef SRC_LANCET_CALLER_GENOTYPER_H_
#define SRC_LANCET_CALLER_GENOTYPER_H_

#include <array>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

extern "C" {
#include "minimap.h"
}

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/caller/variant_support.h"
#include "lancet/hts/cigar_unit.h"
#include "lancet/hts/cigar_utils.h"
#include "lancet/cbdg/read.h"

namespace lancet::caller {

class VariantSet;
class RawVariant;
/*
 * ============================================================================
 * ALIGNMENT PARADIGM SHIFT: MSA vs. Read-to-Haplotype Mapping
 * ============================================================================
 * There is an intentional paradox in the pipeline's alignment parameters:
 * We use highly FORGIVING parameters to build the Contig-Reference MSA, but
 * highly STRICT parameters when mapping raw reads back to the haplotypes.
 * 
 * This is because the expected source of sequence divergence fundamentally
 * shifts between Variant Discovery and Genotyping.
 *
 * 1. THE MSA (Contig vs. Reference) -> Modeling BIOLOGY (The Mapmaker)
 *    - Assumption: The contig is high-confidence. Divergence is true mutation.
 *    - Strategy: FORGIVING. We use cheap gap extensions and high mismatch 
 *      tolerance to force the algorithm to stretch across massive biological 
 *      indels and dense MNVs without soft-clipping.
 *    - Goal: Discover and catalog the variant by keeping the alignment intact.
 *
 * 2. READ-TO-HAPLOTYPE -> Modeling PHYSICS (The GPS)
 *    - Assumption: Raw reads are noisy, but the biological variants are ALREADY 
 *      baked into the assembled haplotype sequences. 
 *    - Divergence: Sequencer error, adapter garbage, or chimeric artifacts.
 *    - Strategy: STRICT. A read should perfectly match its parent haplotype 
 *      with zero biological gaps. We use heavy gap/mismatch penalties to 
 *      punish sequencer noise, force soft-clipping of garbage read-tails, and 
 *      prevent "allele bleeding" (noisy reads mapping to the wrong allele).
 *    - Goal: Accurately segregate read support to calculate clean VAFs.
 * ============================================================================
 */
// ============================================================================
// Genotyper: minimap2-based read-to-haplotype alignment for genotyping
//
// Scoring parameters are custom for Illumina read-to-contig realignment,
// NOT the standard 'sr' preset. See the anonymous namespace block in
//     genotyper.cpp for raw scoring matrix and penalty rationales.
//
// Isolation boundary: AssignReadToAlleles() encapsulates the alignment 
// engine. Everything downstream (AddToTable, VariantSupport) is decoupled.
// ============================================================================
class Genotyper {
 public:
  Genotyper();

  void SetNumSamples(const usize num_samples) { mNumSamples = num_samples; }

  using Reads = absl::Span<const cbdg::Read>;
  using Haplotypes = absl::Span<const std::string>;
  using Result = absl::flat_hash_map<const RawVariant*, SupportArray>;

  [[nodiscard]] auto Genotype(Haplotypes hap_seqs, Reads qry_reads, const VariantSet& variant_set) -> Result;

 private:
  // ──────────────────────────────────────────────────────────────────────────
  // Minimap2 RAII wrappers
  // ──────────────────────────────────────────────────────────────────────────
  struct MmIdxDeleter {
    void operator()(mm_idx_t* idx) noexcept { mm_idx_destroy(idx); }
  };

  struct MmTbufDeleter {
    void operator()(mm_tbuf_t* tbuf) noexcept { mm_tbuf_destroy(tbuf); }
  };

  using MappingOpts = std::unique_ptr<mm_mapopt_t>;
  using IndexingOpts = std::unique_ptr<mm_idxopt_t>;
  using ThreadBuffer = std::unique_ptr<mm_tbuf_t, MmTbufDeleter>;
  using Minimap2Index = std::unique_ptr<mm_idx_t, MmIdxDeleter>;

  static constexpr usize REF_HAP_IDX = 0;

  // ──────────────────────────────────────────────────────────────────────────
  // Outer Class Variables Block (Sorted by descending size: 24B -> 8B -> 4B)
  // ──────────────────────────────────────────────────────────────────────────
  std::vector<Minimap2Index> mIndices;
  std::vector<std::vector<u8>> mEncodedHaplotypes;  // numeric-encoded haplotypes for local scoring
  MappingOpts mMappingOpts = std::make_unique<mm_mapopt_t>();
  IndexingOpts mIndexingOpts = std::make_unique<mm_idxopt_t>();
  ThreadBuffer mThreadBuffer = ThreadBuffer(mm_tbuf_init());
  usize mNumSamples = 0;

  // ──────────────────────────────────────────────────────────────────────────
  // Alignment result from mm_map for a single read-to-haplotype alignment
  // ──────────────────────────────────────────────────────────────────────────
  struct Mm2AlnResult {
    std::vector<hts::CigarUnit> cigar;
    f64 identity = 0.0;     // gap-compressed identity
    usize hap_idx = 0;
    i32 score = 0;          // DP alignment score
    i32 ref_start = 0;      // 0-based start on haplotype (critical for local score offset)
    i32 ref_end = 0;        // 0-based end on haplotype
  };

  // ──────────────────────────────────────────────────────────────────────────
  // Alignment Abstraction Boundary
  //
  // AssignReadToAlleles() completely encapsulates the mapping strategy
  // using minimap2 per haplotype and dynamic local scoring from CIGARs natively.
  //
  // Everything downstream (AddToTable, VariantSupport, VariantCall, VCF)
  // is strictly decoupled from the alignment engine and stays identical.
  // ──────────────────────────────────────────────────────────────────────────
  struct ReadAlleleAssignment {
    // ── Scoring components for allele assignment ──
    //
    // Each read-haplotype pair produces multiple independent signals. The combined
    // score integrates them to assign the read to its best structurally mapped allele:
    //
    //   combined = (global_score - local_raw_score - sc_penalty) + (local_pbq_score * local_identity)
    //
    // Components:
    //
    //   global_score:    mm_map DP score of the full read→haplotype alignment.
    //                    Captures how well the entire read fits this haplotype natively,
    //                    including flanking contexts.
    //
    //   local_raw_score: The absolute raw substitution matrix score of the variant 
    //                    overlap slice. We explicitly SUBTRACT this from global_score 
    //                    to prevent double-counting the variant when we add the PBQ score!
    //
    //   sc_penalty:      Explicit penalization of soft-clipped read tails. Prevents 
    //                    noisy supplementary mappings from artificially inflating 
    //                    their assignment affinities over cleaner native alignments.
    //
    //   local_pbq_score: PBQ-weighted DP score within the variant region only.
    //                    Scales substitution scores by Phred confidence natively 
    //                    (1 - 10^(-PBQ/10)), analogous to GATK's local PairHMM.
    //
    //   local_identity:  Fraction of exact matches natively inside the variant region.
    //                    Acts as a confidence gate on the local PBQ score. A high
    //                    local_score from a noisy alignment (low identity) is
    //                    mathematically discounted, while one from a clean alignment is trusted.
    //
    //   - Why subtract `local_raw_score`?
    //     The `global_score` natively includes the raw unrestricted matrix alignment cost 
    //     spanning the variant sub-region. If we naively appended `local_pbq_score`, we would 
    //     mathematically double-count the locus mapping two variant weights. Subtracting 
    //     `local_raw_score` carves an exact algebraic "hole" out of the global alignment 
    //     path, allowing us to drop the high-fidelity PBQ-weighted score cleanly into 
    //     that specific locus natively.
    //
    //   - Why subtract `sc_penalty`?
    //     Soft-clipped read tails are unaligned sequence garbage. Minimap2 naturally 
    //     exempts them from the primary DP score tracing. By explicitly calculating and 
    //     subtracting a soft-clip penalty, we aggressively crush the combined global 
    //     scores of chimeric supplementary mappings, structurally preventing partially 
    //     matching noise reads from maliciously winning allele assignments against clean 
    //     end-to-end trace alignments.
    //
    //   - Why (local_pbq_score * local_identity)?
    //     This mathematically penalizes "lucky" alignments in low-complexity boundaries.
    //     A noisy chimeric read might accumulate a high DP score (magnitude) simply 
    //     by traversing a chaotic STR locus. By mapping that raw magnitude against 
    //     the CIGAR exact-match fraction (cleanliness), we enforce a strict confidence 
    //     gate structurally. Perfect alignments (identity = 1.0) retain full PBQ weight; 
    //     fragmented, heavily-gapped alignments (e.g. identity < 0.7) are aggressively 
    //     discounted, natively sinking repeat region artifacts and structural chimeras.
    //
    // ── Folded read position ──
    // Folded read position: min(p, 1−p) where p = variant_query_pos / read_length.
    // 0.0 = variant at read edge, 0.5 = variant at read center.
    //
    // Why fold? Artifacts from 3' quality degradation and 5' soft-clip
    // misalignment cluster at BOTH read ends. With raw positions, these
    // bimodal clusters (e.g. positions 5 and 145 in a 150bp read) average
    // to ~75 — indistinguishable from a true centered variant. Folding maps
    // both ends to the same low-value space, converting the bimodal trap
    // into a unidirectional signal: "Are ALT alleles systematically closer
    // to read edges than REF alleles?" Used for RPCD FORMAT field.
    //
    // ── Edit distance (NM) ──
    // Edit distance (NM) of this read against the REF haplotype (hap_idx=0).
    // Mismatches (under M ops, comparing query vs encoded REF) + insertion
    // bases + deletion bases. Soft clips, hard clips, N-skips excluded per
    // SAM spec. Always measured against the reference haplotype regardless
    // of allele assignment so that ASMD = mean(ALT NM) − mean(REF NM)
    // cancels the variant's own contribution and isolates excess noise.
    //
    // ── Representative base quality ──
    // Representative base quality at this variant for this read.
    //
    // For SNVs: the single base quality at the variant position.
    // For indels: the MINIMUM base quality across all read positions spanning
    //             the variant region (weakest-link summary).
    //
    // Why minimum (not mean/median)?
    //   The confidence in observing a complete indel is bounded by the least
    //   confident base in the region. A 10bp deletion where 9 bases are Q30
    //   and 1 base is Q5 should not be treated as high-confidence.
    //
    // How other callers handle this:
    //   - GATK HaplotypeCaller: PairHMM produces a single per-read likelihood
    //     that integrates all base qualities across the alignment. One read
    //     always contributes exactly one observation regardless of variant size.
    //   - bcftools mpileup: uses the minimum base quality in the indel region
    //     as the representative quality for that read.
    //
    // We follow the bcftools convention since we don't have a full PairHMM.
    // The collapse to a single value here ensures that downstream PL and PBQ
    // computations correctly treat each read as one independent observation.

    // ────────────────────────────────────────────────────────────────────────
    // Extracted Member Variables (Sorted by descending alignment size: 8B -> 4B -> 1B)
    // ────────────────────────────────────────────────────────────────────────
    f64 local_score = 0.0;      // (8B) PBQ-weighted DP score within variant region
    f64 local_identity = 0.0;   // (8B) fraction of exact matches in variant region
    f64 folded_read_pos = 0.0;  // (8B) Used for RPCD FORMAT field

    i32 global_score = 0;       // (4B) mm_map DP score of full read→haplotype alignment
    u32 ref_nm = 0;             // (4B) Edit distance (NM) of this read against the REF haplotype
    AlleleIndex allele;         // (4B) Allele enumerator 

    u8 base_qual_at_var = 0;    // (1B) Representative base quality at this variant for this read

    [[nodiscard]] auto CombinedScore() const -> f64 {
      return static_cast<f64>(global_score) + local_score * local_identity;
    }
  };

  void ResetData(Haplotypes hap_seqs);
  
  using PerVariantAssignment = absl::flat_hash_map<const RawVariant*, ReadAlleleAssignment>;
  [[nodiscard]] auto AssignReadToAlleles(const cbdg::Read& qry_read, const VariantSet& variant_set) -> PerVariantAssignment;

  [[nodiscard]] auto AlignToAllHaplotypes(const cbdg::Read& qry_read) -> std::vector<Mm2AlnResult>;
  
  [[nodiscard]] static auto EncodeSequence(std::string_view hap_seq) -> std::vector<u8>;

  // ──────────────────────────────────────────────────────────────────────────
  // AssignReadToAlleles Internal Core Helpers
  // ──────────────────────────────────────────────────────────────────────────
  auto ComputeRefEditDistance(const std::vector<Mm2AlnResult>& alns,
                              absl::Span<const u8> qry_seq_encoded, 
                              usize qry_read_length) const -> u32;

  auto EvaluateAlignment(const Mm2AlnResult& aln, 
                         const RawVariant& variant,
                         absl::Span<const u8> qry_seq_encoded, 
                         absl::Span<const u8> qry_base_quals, 
                         usize qry_read_length,
                         ReadAlleleAssignment& best, 
                         f64& best_combined) const -> bool;

  auto ExtractHapBounds(const RawVariant& variant, usize aln_hap_idx,
                        i32& out_hap_variant_start, i32& out_hap_variant_length, AlleleIndex& out_allele) const -> bool;

  static void AddToTable(Result& out_vars_table, const cbdg::Read& qry_read, const PerVariantAssignment& allele_assignments);
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_GENOTYPER_H_
