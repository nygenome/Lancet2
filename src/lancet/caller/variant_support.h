#ifndef SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
#define SRC_LANCET_CALLER_VARIANT_SUPPORT_H_

#include <array>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "lancet/base/types.h"

namespace lancet::caller {

// ============================================================================
// AlleleIndex: integer-based allele identifier for multi-allelic support
//
// Unlike the previous binary Allele enum {REF, ALT}, AlleleIndex supports
// arbitrary numbers of alleles:
//   0 = REF, 1 = ALT1, 2 = ALT2, ...
//
// This is the foundation for multi-allelic VCF output (comma-separated ALTs)
// and for Phase 2 graph alignment where a read's path through the POA graph
// directly identifies which allele it supports at each variant position.
// ============================================================================
using AlleleIndex = u8;
static constexpr AlleleIndex REF_ALLELE_IDX = 0;

enum class Strand : bool { FWD, REV };

// ============================================================================
// VariantSupport: per-sample allele evidence aggregator
//
// Collects read-level evidence for each allele at a variant site and computes
// aggregate metrics for VCF FORMAT fields. Designed for multi-allelic sites.
//
//  ┌──────────────────┐
//  │  Read Evidence    │  ← one per aligned read at this variant
//  │  {allele, quals}  │
//  └────────┬─────────┘
//           ▼
//  ┌──────────────────┐
//  │  PerAlleleData    │  ← one per allele (REF, ALT1, ALT2...)
//  │  [fwd_bq, rev_bq │     stored in dense vector indexed by AlleleIndex
//  │   map_quals, ...]│
//  └────────┬─────────┘
//           ▼
//  ┌──────────────────┐
//  │ Aggregation       │  → PL, PBQ, RMQ, SB, GQ for VCF FORMAT
//  └──────────────────┘
//
// The PerAlleleData is stored as std::vector<PerAlleleData> indexed directly
// by AlleleIndex. This is efficient because AlleleIndex is a dense, zero-based
// integer (typically 0-3 for most sites, never more than ~8).
// ============================================================================
class VariantSupport {
 public:
  VariantSupport() = default;

  struct ReadEvidence {
    i64 insert_size;        // template length from original alignment
    f64 aln_score;          // normalized alignment score to the assigned haplotype
    f64 folded_read_pos;    // 0.0=read edge, 0.5=read center (for RPRS)
    u32 rname_hash;         // hash of read name (for dedup)
    u32 ref_nm;             // edit distance to REF haplotype (for ASMD)
    AlleleIndex allele;     // which allele this read supports
    Strand strand;          // forward or reverse strand
    u8 base_qual;           // representative PBQ (min across variant region for indels)
    u8 map_qual;            // original mapping quality of the read
    bool is_soft_clipped;   // soft-clip bases >= 6% of read length in original alignment
    bool is_proper_pair;    // properly paired in original alignment
  };

  void AddEvidence(const ReadEvidence& evidence);

  // ── Per-Allele Accessors ──
  [[nodiscard]] auto FwdCount(AlleleIndex idx) const -> usize;
  [[nodiscard]] auto RevCount(AlleleIndex idx) const -> usize;
  [[nodiscard]] auto TotalAlleleCov(AlleleIndex idx) const -> usize;
  [[nodiscard]] auto TotalSampleCov() const noexcept -> usize;
  [[nodiscard]] auto NumAlleles() const noexcept -> usize { return mAlleleData.size(); }

  // Convenience shims matching the old REF/ALT interface (for VariantCall)
  [[nodiscard]] auto TotalRefCov() const -> usize { return TotalAlleleCov(REF_ALLELE_IDX); }
  [[nodiscard]] auto TotalAltCov() const -> usize;

  // ── Aggregate Metrics for VCF FORMAT Fields ──

  // Posterior base quality: Bayesian aggregation of per-read error probabilities
  // across all reads supporting this allele (Edgar & Flyvbjerg 2014).
  //
  // Given N reads with Phred qualities Q_i:
  //   ε_i = 10^(-Q_i/10)
  //   log_err  = Σ log10(ε_i)
  //   log_ok   = Σ log10(1 - ε_i)
  //   posterior = 10^log_err / (10^log_err + 10^log_ok)
  //   PBQ = -10 * log10(posterior), capped at 255
  [[nodiscard]] auto PosteriorBaseQual(AlleleIndex idx) const -> f64;

  // RMS mapping quality: sqrt(mean(mapq_i^2)) for reads supporting this allele.
  // Per-allele FORMAT field (RMQ), follows the standard samtools/GATK convention.
  [[nodiscard]] auto RmsMappingQual(AlleleIndex idx) const -> f64;

  // Phred-scaled Fisher's exact test strand bias: -10*log10(p_value) from a 2×2
  // contingency table of REF/ALT × FWD/REV strand counts. Artifacts concentrate
  // on one strand and produce high scores; true variants score near 0.
  // This is a single per-sample metric (Number=1), not per-allele.
  [[nodiscard]] auto PhredStrandBias() const -> f64;

  // Mean normalized alignment score for reads supporting this allele.
  [[nodiscard]] auto MeanAlnScore(AlleleIndex idx) const -> f64;

  // Soft Clip Asymmetry (SCA FORMAT field): fraction of soft-clipped reads in ALT
  // minus fraction in REF. Detects unresolved larger variant events masquerading
  // as smaller local calls.
  //   SCA = (alt_sc / alt_total) - (ref_sc / ref_total)
  [[nodiscard]] auto SoftClipAsymmetry() const -> f64;

  // Fragment Length Delta (FLD FORMAT field): absolute difference in mean insert
  // sizes between ALT-supporting and REF-supporting read pairs. Large FLD
  // indicates chimeric library artifacts or somatic cfDNA fragment length shifts.
  //   FLD = |mean_alt_isize - mean_ref_isize|
  [[nodiscard]] auto FragLengthDelta() const -> f64;

  // Mapping Quality Rank Sum Z-score (MQRS FORMAT field): Mann-Whitney U test
  // comparing mapping qualities of REF-supporting vs ALT-supporting reads.
  // Uses the original alignment MAPQ (not re-alignment quality).
  //   Negative Z → ALT reads are systematically lower MAPQ (paralogous mismapping)
  //   Positive Z → ALT reads have higher MAPQ (unusual, may indicate REF bias)
  //   Near zero  → no systematic difference (expected for true variants)
  //   100.0      → test not applicable (empty REF or ALT group, or identical MAPQ)
  // GATK uses Z < -12.5 as a hard filter threshold for MQRankSum.
  [[nodiscard]] auto MappingQualRankSumZ() const -> f64;

  // Read Position Rank Sum Z-score (RPRS FORMAT field): Mann-Whitney U test
  // on folded read positions (REF vs ALT). Folded position maps both 5' and
  // 3' read edges to 0.0, centers to 0.5. True variants are uniformly
  // distributed; artifacts from quality degradation or misalignment cluster
  // at read edges, producing a negative Z-score for the ALT group.
  //   Negative Z → ALT allele systematically at read edges (artifact signal)
  //   Near zero  → uniform distribution (expected for true variants)
  //   100.0      → untestable (empty group or identical distributions)
  [[nodiscard]] auto ReadPosRankSumZ() const -> f64;

  // Base Quality Rank Sum Z-score (BQRS FORMAT field): Mann-Whitney U test
  // on per-allele base qualities (REF vs ALT). Detects 8-oxoguanine oxidation
  // artifacts and other chemistry-driven errors where the miscalled base has
  // characteristically low Phred confidence.
  //   Negative Z → ALT allele has systematically lower base quality
  //   100.0      → untestable
  [[nodiscard]] auto BaseQualRankSumZ() const -> f64;

  // Allele-Specific Mismatch Delta (ASMD FORMAT field):
  // mean(ALT NM) − mean(REF NM), where NM is edit distance to the REF haplotype.
  // True variants contribute equally to both groups' NM (the variant itself is
  // an edit against the reference for all reads). Chimeric or paralogously
  // mismapped reads produce high ALT NM (excess random mismatches) while REF
  // reads stay clean, yielding ASMD > 0.
  //   0.0 if either group is empty.
  [[nodiscard]] auto AlleleMismatchDelta() const -> f64;

  // ── Multi-Allelic Genotype Likelihoods ──

  // Computes Phred-scaled genotype likelihoods (PLs) for all possible diploid
  // genotypes at a multi-allelic site.
  //
  // For k alleles (0=REF, 1..k-1=ALTs), there are k*(k+1)/2 diploid genotypes.
  // Uses the standard GATK per-read likelihood model:
  //
  //   P(read | GT=(a1,a2)) = 0.5 * P(read|a1) + 0.5 * P(read|a2)
  //   P(read | allele_a) = (1-ε) if read matches allele a, else ε/(k-1)
  //
  // Genotypes are ordered per VCF 4.3 spec §1.6.2:
  //   GL_index(i,j) = j*(j+1)/2 + i  where i ≤ j
  //
  // Returns: vector of k*(k+1)/2 normalized Phred-scaled likelihoods (min PL=0)
  //
  // References:
  //   - Li, H. (2011). "A statistical framework for SNP calling..."
  //   - Poplin, R. et al. (2018). GATK HaplotypeCaller
  //   - VCF 4.3 specification, Section 1.6.2 (GL index formula)
  [[nodiscard]] auto ComputePLs() const -> std::vector<int>;

  // Genotype Quality (GQ): confidence in the called genotype.
  // Defined as the difference between the second-lowest and lowest PL values,
  // capped at 99. Standard GATK convention.
  //   GQ = second_min(PLs) - min(PLs)    [min is always 0 after normalization]
  // See: https://gatk.broadinstitute.org/hc/en-us/articles/360035531692
  [[nodiscard]] static auto ComputeGQ(const std::vector<int>& pls) -> int;

  // Copy allele data from `src` allele `src_allele` into `dst_allele` slot
  // in this object. Used for multi-allelic merging: each bi-allelic variant
  // has alleles {0=REF, 1=ALT}. When merging N variants at a locus, we remap
  // variant[i]'s ALT(1) → merged allele (i+1).
  void MergeAlleleFrom(const VariantSupport& src, AlleleIndex src_allele, AlleleIndex dst_allele);


 private:
  struct PerAlleleData {
    // Deduplicate reads: a read can support an allele only once per strand.
    // Key = read name hash, value = strand seen.
    absl::flat_hash_map<u32, Strand> name_hashes;

    // Per-read representative base quality at the variant position, split by
    // strand. For indels, this is the MINIMUM PBQ across the variant region
    // (one entry per read, NOT per base position — see AddEvidence comment).
    std::vector<u8> fwd_base_quals;
    std::vector<u8> rev_base_quals;

    // Per-read mapping quality (for RMS RMQ computation).
    std::vector<u8> map_quals;

    // Per-read normalized alignment score (for filtering / annotation).
    std::vector<f64> aln_scores;

    // Count of soft-clipped reads supporting this allele (for SCA FORMAT tag).
    usize soft_clip_count = 0;

    // Insert sizes from properly-paired reads (for FLD FORMAT tag).
    // Only non-zero insert sizes from proper pairs are tracked.
    std::vector<f64> proper_pair_isizes;

    // Folded read positions: min(p, 1−p) for RPRS rank-sum test.
    std::vector<f64> folded_read_positions;

    // Edit distances (NM) against REF haplotype for ASMD delta.
    // Stored as f64 for mean computation.
    std::vector<f64> ref_nm_values;
  };

  // Dense vector indexed by AlleleIndex: mAlleleData[0]=REF, [1]=ALT1, ...
  std::vector<PerAlleleData> mAlleleData;

  // Grow the vector to accommodate a new allele index.
  void EnsureAlleleSlot(AlleleIndex idx);
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
