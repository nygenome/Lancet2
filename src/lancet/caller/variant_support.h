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

  // Evidence from a single read's alignment at this variant
  struct ReadEvidence {
    u32 rname_hash;       // hash of read name (for dedup)
    AlleleIndex allele;   // which allele this read supports
    Strand strand;        // forward or reverse strand
    u8 base_qual;         // representative PBQ (min across variant region for indels)
    u8 map_qual;          // original mapping quality of the read
    f64 aln_score;        // normalized alignment score to the assigned haplotype
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
  // Follows the standard samtools/GATK INFO/RMQ convention.
  [[nodiscard]] auto RmsMappingQual(AlleleIndex idx) const -> f64;

  // Strand bias: fraction of forward-strand reads (fwd / total) for this allele.
  [[nodiscard]] auto StrandBiasRatio(AlleleIndex idx) const -> f64;

  // Mean normalized alignment score for reads supporting this allele.
  [[nodiscard]] auto MeanAlnScore(AlleleIndex idx) const -> f64;

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
  };

  // Dense vector indexed by AlleleIndex: mAlleleData[0]=REF, [1]=ALT1, ...
  std::vector<PerAlleleData> mAlleleData;

  // Grow the vector to accommodate a new allele index.
  void EnsureAlleleSlot(AlleleIndex idx);
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
