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

namespace genotyper_detail {

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

}  // namespace genotyper_detail

// ============================================================================
// Genotyper: minimap2-based read-to-haplotype alignment for genotyping
//
// Scoring parameters are custom for Illumina read-to-contig realignment,
// NOT the standard 'sr' preset. See genotyper_detail for rationale.
//
// Phase 2 boundary: AssignReadToAlleles() will be replaced by graph DP
// (gssw/vg). Everything downstream (AddToTable, VariantSupport) stays.
// ============================================================================
class Genotyper {
 public:
  Genotyper();

  void SetNumSamples(const usize num_samples) { mNumSamples = num_samples; }

  using Reads = absl::Span<const cbdg::Read>;
  using Haplotypes = absl::Span<const std::string>;
  using Result = absl::flat_hash_map<const RawVariant*, SupportArray>;

  [[nodiscard]] auto Genotype(Haplotypes haps, Reads reads, const VariantSet& vset) -> Result;

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

  usize mNumSamples = 0;
  std::vector<Minimap2Index> mIndices;
  std::vector<std::vector<u8>> mEncodedHaplotypes;  // numeric-encoded haplotypes for local scoring
  MappingOpts mMappingOpts = std::make_unique<mm_mapopt_t>();
  IndexingOpts mIndexingOpts = std::make_unique<mm_idxopt_t>();
  ThreadBuffer mThreadBuffer = ThreadBuffer(mm_tbuf_init());

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
    bool is_full_match = false;  // read aligns end-to-end with high identity
  };

  // ──────────────────────────────────────────────────────────────────────────
  // Phase 2 Abstraction Boundary
  //
  // AssignReadToAlleles() encapsulates the alignment strategy:
  //   Phase 1: minimap2 mm_map per haplotype → local scoring from CIGAR
  //   Phase 2: replace with graph DP against POA graph (gssw/vg)
  //
  // Everything downstream (AddToTable, VariantSupport, VariantCall, VCF)
  // is decoupled from alignment strategy and stays identical.
  // ──────────────────────────────────────────────────────────────────────────
  struct ReadAlleleAssignment {
    // ── Scoring components for allele assignment ──
    //
    // Each read-haplotype pair produces three independent signals. The combined
    // score integrates all three to assign the read to its best allele:
    //
    //   combined = global_score + local_score * local_identity
    //
    // Three-signal model:
    //
    //   global_score:    mm_map DP score of the full read→haplotype alignment.
    //                    Captures how well the entire read fits this haplotype,
    //                    including flanking context. Prevents reads with isolated
    //                    matches from incorrectly winning allele assignment.
    //
    //   local_score:     PBQ-weighted DP score within the variant region only.
    //                    Each base's substitution score is scaled by its Phred
    //                    confidence (1 - 10^(-PBQ/10)), so low-quality bases
    //                    contribute less. This is analogous to GATK's PairHMM
    //                    which integrates PBQ into per-read log-likelihoods.
    //
    //   local_identity:  Fraction of exact matches in the variant region.
    //                    Acts as a confidence gate on the local score. A high
    //                    local_score from a noisy alignment (low identity) is
    //                    discounted, while one from a clean alignment is trusted.
    //
    // Why this formula instead of ad-hoc weights?
    //
    //   - GATK HaplotypeCaller: PairHMM produces a single integrated
    //     log-likelihood P(read | haplotype). We don't have a PairHMM, but
    //     global_score approximates log P(read | haplotype) in score space.
    //
    //   - Strelka2: enumerates candidate alignments and takes the maximum
    //     alignment-specific likelihood. Our per-haplotype mm_map is analogous.
    //
    //   - The multiplicative local_score * local_identity term models the
    //     quality-of-evidence at the variant: local_score measures magnitude,
    //     local_identity measures cleanliness. Their product is a
    //     quality-weighted local contribution — if identity is 1.0 (perfect),
    //     local_score contributes fully; if identity is 0.5 (noisy), it's
    //     halved. This naturally handles repeat regions and alignment artifacts
    //     where high DP scores coexist with poor identity.
    //
    //   - Phase 2 (graph DP) will replace this with a single graph-alignment
    //     log-likelihood that naturally integrates all three signals.
    f64 local_score;        // PBQ-weighted DP score within variant region
    f64 local_identity;     // fraction of exact matches in variant region

    /// Folded read position: min(p, 1−p) where p = variant_query_pos / read_length.
    /// 0.0 = variant at read edge, 0.5 = variant at read center.
    ///
    /// Why fold? Artifacts from 3' quality degradation and 5' soft-clip
    /// misalignment cluster at BOTH read ends. With raw positions, these
    /// bimodal clusters (e.g. positions 5 and 145 in a 150bp read) average
    /// to ~75 — indistinguishable from a true centered variant. Folding maps
    /// both ends to the same low-value space, converting the bimodal trap
    /// into a unidirectional signal: "Are ALT alleles systematically closer
    /// to read edges than REF alleles?" Used for RPCD FORMAT field.
    f64 folded_read_pos = 0.0;

    i32 global_score;       // mm_map DP score of full read→haplotype alignment

    /// Edit distance (NM) of this read against the REF haplotype (hap_idx=0).
    /// Mismatches (under M ops, comparing query vs encoded REF) + insertion
    /// bases + deletion bases. Soft clips, hard clips, N-skips excluded per
    /// SAM spec. Always measured against the reference haplotype regardless
    /// of allele assignment so that ASMD = mean(ALT NM) − mean(REF NM)
    /// cancels the variant's own contribution and isolates excess noise.
    u32 ref_nm = 0;

    AlleleIndex allele;

    [[nodiscard]] auto CombinedScore() const -> f64 {
      return static_cast<f64>(global_score) + local_score * local_identity;
    }

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
    u8 base_qual_at_var;
  };

  using PerVariantAssignment = absl::flat_hash_map<const RawVariant*, ReadAlleleAssignment>;
  [[nodiscard]] auto AssignReadToAlleles(const cbdg::Read& read,
                                          const VariantSet& vset) -> PerVariantAssignment;

  [[nodiscard]] auto AlignToAllHaplotypes(const cbdg::Read& read) -> std::vector<Mm2AlnResult>;

  void ResetData(Haplotypes seqs);

  static void AddToTable(Result& rslt, const cbdg::Read& read,
                          const PerVariantAssignment& assignments);

  [[nodiscard]] static auto EncodeSequence(std::string_view seq) -> std::vector<u8>;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_GENOTYPER_H_
