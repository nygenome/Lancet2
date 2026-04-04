#include "lancet/caller/variant_support.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "lancet/base/types.h"
#include "lancet/hts/phred_quality.h"

namespace lancet::caller {

// ============================================================================
// AddEvidence
// ============================================================================
void VariantSupport::AddEvidence(const ReadEvidence& evidence) {
  EnsureAlleleSlot(evidence.allele);
  auto& data = mAlleleData[evidence.allele];

  // Deduplicate by read name hash: only first-seen strand counts.
  const auto [iter, inserted] = data.name_hashes.try_emplace(evidence.rname_hash, evidence.strand);
  if (!inserted) {
    return;
  }

  // Each read contributes exactly one representative base quality value.
  // For indels, this is already collapsed to min(PBQ) across the variant region
  // by the genotyper (see ReadAlleleAssignment::base_qual_at_var comment).
  if (evidence.strand == Strand::FWD) {
    data.fwd_base_quals.push_back(evidence.base_qual);
  } else {
    data.rev_base_quals.push_back(evidence.base_qual);
  }

  data.map_quals.push_back(evidence.map_qual);
  data.aln_scores.push_back(evidence.aln_score);
}

// ============================================================================
// MergeAlleleFrom — copy one allele's read data from src into dst_allele slot.
//
// Used by the multi-allelic VariantCall constructor to build a unified
// VariantSupport from N bi-allelic per-variant supports:
//
//   variant[0]: {REF=allele0, ALT=allele1}  → merge REF→0, ALT→1
//   variant[1]: {REF=allele0, ALT=allele1}  → merge ALT→2
//   variant[2]: {REF=allele0, ALT=allele1}  → merge ALT→3
//
// Deduplicates by read name hash (same read can't count twice).
// ============================================================================
void VariantSupport::MergeAlleleFrom(const VariantSupport& src, const AlleleIndex src_allele,
                                     const AlleleIndex dst_allele) {
  if (src_allele >= src.mAlleleData.size()) return;
  EnsureAlleleSlot(dst_allele);

  const auto& src_data = src.mAlleleData[src_allele];
  auto& dst_data = mAlleleData[dst_allele];

  for (const auto& [rname_hash, strand] : src_data.name_hashes) {
    const auto [iter, inserted] = dst_data.name_hashes.try_emplace(rname_hash, strand);
    if (!inserted) continue;  // already seen this read in dst
  }

  // Append PBQ, RMQ, and alignment scores (reads already deduped above by name_hashes,
  // but the qual vectors are append-only during AddEvidence, so we append here too).
  dst_data.fwd_base_quals.insert(dst_data.fwd_base_quals.end(),
                                 src_data.fwd_base_quals.begin(), src_data.fwd_base_quals.end());
  dst_data.rev_base_quals.insert(dst_data.rev_base_quals.end(),
                                 src_data.rev_base_quals.begin(), src_data.rev_base_quals.end());
  dst_data.map_quals.insert(dst_data.map_quals.end(),
                            src_data.map_quals.begin(), src_data.map_quals.end());
  dst_data.aln_scores.insert(dst_data.aln_scores.end(),
                             src_data.aln_scores.begin(), src_data.aln_scores.end());
}

// ============================================================================
// Per-Allele Accessors
// ============================================================================
auto VariantSupport::FwdCount(const AlleleIndex idx) const -> usize {
  if (idx >= mAlleleData.size()) {
    return 0;
  }
  return mAlleleData[idx].fwd_base_quals.size();
}

auto VariantSupport::RevCount(const AlleleIndex idx) const -> usize {
  if (idx >= mAlleleData.size()) {
    return 0;
  }
  return mAlleleData[idx].rev_base_quals.size();
}

auto VariantSupport::TotalAlleleCov(const AlleleIndex idx) const -> usize {
  return FwdCount(idx) + RevCount(idx);
}

auto VariantSupport::TotalSampleCov() const noexcept -> usize {
  usize total = 0;
  for (usize i = 0; i < mAlleleData.size(); ++i) {
    total += mAlleleData[i].fwd_base_quals.size() + mAlleleData[i].rev_base_quals.size();
  }
  return total;
}

auto VariantSupport::TotalAltCov() const -> usize {
  usize total = 0;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    total += TotalAlleleCov(static_cast<AlleleIndex>(i));
  }
  return total;
}

// ============================================================================
// PosteriorBaseQual (PBQ FORMAT field)
//
// Edgar & Flyvbjerg (2014): Bayesian aggregation of per-read error probs.
//
// Given N reads supporting allele `a`, each with Phred quality Q_i:
//
//   ε_i = 10^(-Q_i / 10)
//
//   log10(P_err_all)  = Σ log10(ε_i)           // all reads are wrong
//   log10(P_ok_all)   = Σ log10(1 - ε_i)       // all reads are correct
//
//   posterior_err = P_err_all / (P_err_all + P_ok_all)
//
// In log10 space (for numerical stability via log-sum-exp):
//   log10(posterior_err) = log_err - log10(10^log_err + 10^log_ok)
//
//   PBQ = -10 * log10(posterior_err), capped at 255
//
// Intuition: if two reads agree with PBQ=30 each, the combined posterior
// base quality is much higher than 30, reflecting increased confidence.
// ============================================================================
auto VariantSupport::PosteriorBaseQual(const AlleleIndex idx) const -> f64 {
  if (idx >= mAlleleData.size()) {
    return 0.0;
  }

  const auto& data = mAlleleData[idx];
  // Per-read representative base quality at the variant position, split by
  // strand. For indels, this is the MINIMUM PBQ across the variant region
  // (one entry per read, NOT per base position).
  const auto& fwd_bq = data.fwd_base_quals;
  const auto& rev_bq = data.rev_base_quals;

  if (fwd_bq.empty() && rev_bq.empty()) {
    return 0.0;
  }

  f64 log_err = 0.0;  // Σ log10(ε_i)
  f64 log_ok = 0.0;   // Σ log10(1 - ε_i)

  const auto accumulate_bq = [&log_err, &log_ok](const std::vector<u8>& quals) {
    for (const auto qual : quals) {
      const f64 eps = hts::PhredToErrorProb(qual);
      log_err += std::log10(std::max(eps, 1e-300));
      log_ok += std::log10(std::max(1.0 - eps, 1e-300));
    }
  };

  accumulate_bq(fwd_bq);
  accumulate_bq(rev_bq);

  // log-sum-exp: log10(10^a + 10^b) = max(a,b) + log10(1 + 10^(min-max))
  const f64 max_log = std::max(log_err, log_ok);
  const f64 log_sum = max_log + std::log10(1.0 + std::pow(10.0, std::min(log_err, log_ok) - max_log));
  const f64 log_posterior_err = log_err - log_sum;
  const f64 posterior_bq = -10.0 * log_posterior_err;

  static constexpr f64 MAX_BQ = 255.0;
  return std::min(posterior_bq, MAX_BQ);
}

// ============================================================================
// RmsMappingQual (RMQ FORMAT field)
//
// Root-mean-square of the mapping qualities of all reads assigned to this
// allele. This follows the standard INFO/RMQ convention (samtools, GATK).
//
//   RMQ = sqrt( Σ(mapq_i²) / N )
// ============================================================================
auto VariantSupport::RmsMappingQual(const AlleleIndex idx) const -> f64 {
  if (idx >= mAlleleData.size() || mAlleleData[idx].map_quals.empty()) {
    return 0.0;
  }

  const auto& mqs = mAlleleData[idx].map_quals;
  const f64 sum_sq = std::accumulate(mqs.cbegin(), mqs.cend(), 0.0, [](const f64 acc, const u8 mq) {
    return acc + static_cast<f64>(mq) * static_cast<f64>(mq);
  });

  return std::sqrt(sum_sq / static_cast<f64>(mqs.size()));
}

// ============================================================================
// StrandBiasRatio (SB FORMAT field)
//
// Simple strand bias metric: fraction of forward-strand reads for this allele.
// A value of 0.5 indicates no bias. Values near 0 or 1 indicate strong bias.
//
//   SB = fwd_count / total_count
// ============================================================================
auto VariantSupport::StrandBiasRatio(const AlleleIndex idx) const -> f64 {
  const auto total = TotalAlleleCov(idx);
  if (total == 0) {
    return 0.0;
  }
  return static_cast<f64>(FwdCount(idx)) / static_cast<f64>(total);
}

// ============================================================================
// MeanAlnScore
// ============================================================================
auto VariantSupport::MeanAlnScore(const AlleleIndex idx) const -> f64 {
  if (idx >= mAlleleData.size() || mAlleleData[idx].aln_scores.empty()) {
    return 0.0;
  }

  const auto& scores = mAlleleData[idx].aln_scores;
  const f64 sum = std::accumulate(scores.cbegin(), scores.cend(), 0.0);
  return sum / static_cast<f64>(scores.size());
}

// ============================================================================
// ComputePLs — Multi-Allelic Phred-Scaled Genotype Likelihoods
//
// Standard GATK per-read genotype likelihood model for diploid organisms.
//
// For k alleles and a read supporting allele `s` with error probability ε:
//
//   P(read=s | true_allele=a) = (1-ε)       if s == a   (read matches)
//                              = ε/(k-1)     if s != a   (error)
//
//   P(read | GT=(a1,a2)) = 0.5 * P(read|a1) + 0.5 * P(read|a2)
//
//   P(Data | GT) = ∏_reads P(read | GT)     (independence)
//
// Genotype ordering follows VCF 4.3 §1.6.2:
//   GL_index(i,j) = j*(j+1)/2 + i   where i ≤ j
//
//   Bi-allelic (k=2):  {0/0, 0/1, 1/1}
//   Tri-allelic (k=3): {0/0, 0/1, 1/1, 0/2, 1/2, 2/2}
//
// References:
//   Li, H. (2011). "A statistical framework for SNP calling..."
//   Poplin, R. et al. (2018). GATK HaplotypeCaller.
// ============================================================================

namespace {

/// For a single read supporting `read_allele` with base quality `base_qual`,
/// compute P(read | true_allele = a) for every allele a in [0, num_alleles).
///
/// Returns a vector of conditional probabilities, one per allele:
///   result[a] = (1 - ε)      if a == read_allele  (read matches true allele)
///             = ε / (k - 1)  if a != read_allele  (read is a sequencing error)
///
/// This is precomputed once per read [O(k)], then reused across all
/// O(k²/2) genotype pairs — avoiding redundant recomputation.
auto ReadLikelihoodsPerAllele(const int read_allele, const u8 base_qual,
                              const int num_alleles) -> std::vector<f64> {
  const f64 error_prob = hts::PhredToErrorProb(base_qual);
  const f64 match_prob = 1.0 - error_prob;
  const f64 mismatch_prob = error_prob / std::max(1, num_alleles - 1);

  std::vector<f64> likelihoods(num_alleles, mismatch_prob);
  likelihoods[read_allele] = match_prob;
  return likelihoods;
}

/// Accumulate one read's contribution into the genotype log-likelihoods.
///
/// For each genotype (allele_a, allele_b) where a ≤ b:
///   P(read | GT=(a,b)) = 0.5 * P(read|a) + 0.5 * P(read|b)
///   log_GL[gt] += log10(P(read | GT))
///
/// Complexity: O(k²/2) per read, with per-allele likelihoods already computed.
void AccumulateReadIntoGLs(const std::vector<f64>& read_lks,
                           const int num_alleles,
                           std::vector<f64>& gt_log_likelihoods) {
  for (int allele_b = 0; allele_b < num_alleles; ++allele_b) {
    for (int allele_a = 0; allele_a <= allele_b; ++allele_a) {
      const int gt_idx = allele_b * (allele_b + 1) / 2 + allele_a;

      // Diploid mixture: equal probability of inheriting from either chromosome
      const f64 diploid_prob = 0.5 * read_lks[allele_a] + 0.5 * read_lks[allele_b];

      // Guard against log10(0); floor at 10^-300 ≈ -300 in log space
      gt_log_likelihoods[gt_idx] += std::log10(std::max(diploid_prob, 1e-300));
    }
  }
}

/// Convert raw log-likelihoods to normalized, capped Phred-scaled PLs.
/// The best genotype gets PL=0; others are scaled relative to it.
auto NormalizeLogLikelihoodsToPLs(const std::vector<f64>& gt_log_lks) -> std::vector<int> {
  static constexpr f64 PL_CAP = static_cast<f64>(std::numeric_limits<int>::max() / 2);

  const f64 best_ll = *std::max_element(gt_log_lks.cbegin(), gt_log_lks.cend());

  std::vector<int> pls(gt_log_lks.size());
  for (usize idx = 0; idx < gt_log_lks.size(); ++idx) {
    const f64 raw_pl = -10.0 * (gt_log_lks[idx] - best_ll);
    pls[idx] = static_cast<int>(std::round(std::min(raw_pl, PL_CAP)));
  }
  return pls;
}

}  // namespace

auto VariantSupport::ComputePLs() const -> std::vector<int> {
  const auto num_alleles = static_cast<int>(mAlleleData.size());
  if (num_alleles == 0) return {};

  const auto num_genotypes = num_alleles * (num_alleles + 1) / 2;
  std::vector<f64> gt_log_lks(num_genotypes, 0.0);

  // Walk every read once. Fwd and rev strands use identical math —
  // strand information is already captured in the SB (strand bias) field.
  for (int allele_idx = 0; allele_idx < num_alleles; ++allele_idx) {
    const auto& allele_data = mAlleleData[allele_idx];

    for (const auto base_qual : allele_data.fwd_base_quals) {
      const auto read_lks = ReadLikelihoodsPerAllele(allele_idx, base_qual, num_alleles);
      AccumulateReadIntoGLs(read_lks, num_alleles, gt_log_lks);
    }

    for (const auto base_qual : allele_data.rev_base_quals) {
      const auto read_lks = ReadLikelihoodsPerAllele(allele_idx, base_qual, num_alleles);
      AccumulateReadIntoGLs(read_lks, num_alleles, gt_log_lks);
    }
  }

  return NormalizeLogLikelihoodsToPLs(gt_log_lks);
}

// ============================================================================
// ComputeGQ — Genotype Quality
//
// Standard GATK convention: GQ is the second-smallest PL value.
// After normalization, the smallest PL is always 0 (the called genotype).
// GQ = second-smallest PL, capped at 99.
//
// See: https://gatk.broadinstitute.org/hc/en-us/articles/360035531692
// ============================================================================
auto VariantSupport::ComputeGQ(const std::vector<int>& pls) -> int {
  if (pls.size() < 2) {
    return 0;
  }

  // Find minimum and second-minimum PL values
  int min1 = std::numeric_limits<int>::max();
  int min2 = std::numeric_limits<int>::max();
  for (const int val : pls) {
    if (val < min1) {
      min2 = min1;
      min1 = val;
    } else if (val < min2) {
      min2 = val;
    }
  }

  // After normalization, min1 should be 0. GQ = min2 - min1 = min2.
  static constexpr int MAX_GQ = 99;
  return std::min(min2 - min1, MAX_GQ);
}

// ============================================================================
// EnsureAlleleSlot
// ============================================================================
void VariantSupport::EnsureAlleleSlot(const AlleleIndex idx) {
  if (idx >= mAlleleData.size()) {
    mAlleleData.resize(idx + 1);
  }
}

}  // namespace lancet::caller
