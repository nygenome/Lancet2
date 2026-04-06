#include "lancet/caller/variant_support.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "lancet/base/mann_whitney.h"
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

  // Track soft-clip status from original alignment (for SCA FORMAT tag)
  if (evidence.is_soft_clipped) {
    ++data.soft_clip_count;
  }

  // Track insert sizes from properly-paired reads (for FLD FORMAT tag)
  if (evidence.is_proper_pair && evidence.insert_size != 0) {
    data.proper_pair_isizes.push_back(std::abs(static_cast<f64>(evidence.insert_size)));
  }

  // Track folded read position (for RPCD FORMAT tag)
  data.folded_read_positions.push_back(evidence.folded_read_pos);

  // Track edit distance to REF haplotype (for ASMD FORMAT tag)
  data.ref_nm_values.push_back(static_cast<f64>(evidence.ref_nm));
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
  dst_data.soft_clip_count += src_data.soft_clip_count;
  dst_data.proper_pair_isizes.insert(dst_data.proper_pair_isizes.end(),
                                     src_data.proper_pair_isizes.begin(), src_data.proper_pair_isizes.end());
  dst_data.folded_read_positions.insert(dst_data.folded_read_positions.end(),
                                        src_data.folded_read_positions.begin(), src_data.folded_read_positions.end());
  dst_data.ref_nm_values.insert(dst_data.ref_nm_values.end(),
                                src_data.ref_nm_values.begin(), src_data.ref_nm_values.end());
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
// RawPosteriorBaseQual (used to compute NPBQ FORMAT field)
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
//   raw PBQ = -10 * log10(posterior_err)
//
// Returns the raw uncapped value. The caller (BuildFormatFields) divides
// by allele_depth to produce NPBQ (per-read quality) for VCF output.
// ============================================================================
auto VariantSupport::RawPosteriorBaseQual(const AlleleIndex idx) const -> f64 {
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

  return posterior_bq;
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
// StrandBiasLogOR (SB FORMAT field)
//
// Natural log odds ratio for strand bias with Haldane correction (+1 to
// all cells). Coverage-invariant: measures the effect size of strand
// imbalance, not statistical significance.
//
//        |  FWD  |  REV  |
//   -----|-------|-------|
//   REF  | rf    | rr    |
//   ALT  | af    | ar    |      (ALT = sum across all non-REF alleles)
//
//   SB = ln( ((rf+1)(ar+1)) / ((rr+1)(af+1)) )
//
// Positive: ALT enriched on forward strand relative to REF.
// Negative: ALT enriched on reverse strand relative to REF.
// Near 0.0: balanced strands (no bias).
//
// Coverage stability: a 3:1 strand ratio produces SB ≈ 1.1 at any depth
// from 20× to 2000× (vs Fisher SB ranging 3.5→290 over the same range).
// ============================================================================
auto VariantSupport::StrandBiasLogOR() const -> f64 {
  // REF strand counts
  const auto ref_fwd = static_cast<int>(FwdCount(REF_ALLELE_IDX));
  const auto ref_rev = static_cast<int>(RevCount(REF_ALLELE_IDX));

  // ALT strand counts (summed across all non-REF alleles)
  int alt_fwd = 0;
  int alt_rev = 0;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    alt_fwd += static_cast<int>(mAlleleData[i].fwd_base_quals.size());
    alt_rev += static_cast<int>(mAlleleData[i].rev_base_quals.size());
  }

  // Haldane correction: +1 to all cells handles zero-count edge cases
  // without sentinels. When all counts are zero, returns 0.0 (no bias).
  const auto rf1 = static_cast<f64>(ref_fwd + 1);
  const auto rr1 = static_cast<f64>(ref_rev + 1);
  const auto af1 = static_cast<f64>(alt_fwd + 1);
  const auto ar1 = static_cast<f64>(alt_rev + 1);
  return std::log((rf1 * ar1) / (rr1 * af1));
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
// SoftClipAsymmetry (SCA FORMAT field)
//
// Fraction of soft-clipped reads in ALT minus fraction in REF.
//   SCA = (alt_sc/alt_total) - (ref_sc/ref_total)
//
// Positive values indicate ALT reads are more often soft-clipped than REF,
// which may flag unresolved larger structural events masquerading as smaller
// local variants. Soft-clip status comes from the original whole-genome
// alignment CIGAR (>= 6% of read length), not the genotyper re-alignment.
// ============================================================================
auto VariantSupport::SoftClipAsymmetry() const -> f64 {
  // ALT soft-clip fraction (summed across all non-REF alleles)
  usize alt_sc = 0;
  usize alt_total = 0;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    alt_sc += mAlleleData[i].soft_clip_count;
    alt_total += TotalAlleleCov(static_cast<AlleleIndex>(i));
  }

  // REF soft-clip fraction
  const auto ref_sc = (REF_ALLELE_IDX < mAlleleData.size()) ? mAlleleData[REF_ALLELE_IDX].soft_clip_count : usize{0};
  const auto ref_total = TotalAlleleCov(REF_ALLELE_IDX);

  const f64 alt_frac = alt_total > 0 ? static_cast<f64>(alt_sc) / static_cast<f64>(alt_total) : 0.0;
  const f64 ref_frac = ref_total > 0 ? static_cast<f64>(ref_sc) / static_cast<f64>(ref_total) : 0.0;
  return alt_frac - ref_frac;
}

// ============================================================================
// FragLengthDelta (FLD FORMAT field)
//
// Absolute difference in mean insert sizes between ALT-supporting and
// REF-supporting properly-paired reads.
//   FLD = |mean_alt_isize - mean_ref_isize|
//
// Large FLD indicates chimeric library artifacts (artificial bridging) or
// somatic cfDNA fragment length shifts. Insert sizes come from the original
// whole-genome alignment (bam1_t::core.isize), not the genotyper re-alignment.
// Only non-zero insert sizes from properly-paired reads are included.
// ============================================================================
auto VariantSupport::FragLengthDelta() const -> f64 {
  // ALT mean insert size (summed across all non-REF alleles)
  f64 alt_isize_sum = 0.0;
  usize alt_pairs = 0;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    for (const auto isize : mAlleleData[i].proper_pair_isizes) {
      alt_isize_sum += isize;
      ++alt_pairs;
    }
  }

  // REF mean insert size
  f64 ref_isize_sum = 0.0;
  usize ref_pairs = 0;
  if (REF_ALLELE_IDX < mAlleleData.size()) {
    for (const auto isize : mAlleleData[REF_ALLELE_IDX].proper_pair_isizes) {
      ref_isize_sum += isize;
      ++ref_pairs;
    }
  }

  const f64 alt_mean = alt_pairs > 0 ? alt_isize_sum / static_cast<f64>(alt_pairs) : 0.0;
  const f64 ref_mean = ref_pairs > 0 ? ref_isize_sum / static_cast<f64>(ref_pairs) : 0.0;
  return std::abs(alt_mean - ref_mean);
}

// ============================================================================
// MappingQualCohenD (MQCD FORMAT field)
//
// Coverage-normalized effect size comparing mapping qualities of REF vs
// ALT reads. Uses Mann-Whitney U Z/√N to remove √N power amplification.
//
// A mild ALT MAPQ depression (2 units lower) produces MQCD ≈ −0.34 at any
// depth from 20× to 2000× (vs raw Z-score ranging −1.5 to −14.9).
//
// Pools all ALT alleles into a single group (REF vs all-ALT) because the
// test is about whether ALT reads as a class are mismapped.
//
// Returns 0.0 if either group is empty or all MAPQs are identical.
// ============================================================================
auto VariantSupport::MappingQualCohenD() const -> f64 {
  // Collect REF mapping qualities
  absl::Span<const u8> ref_mqs;
  if (REF_ALLELE_IDX < mAlleleData.size()) {
    ref_mqs = absl::MakeConstSpan(mAlleleData[REF_ALLELE_IDX].map_quals);
  }

  // Pool all ALT mapping qualities into a single group
  std::vector<u8> alt_mqs;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    const auto& mqs = mAlleleData[i].map_quals;
    alt_mqs.insert(alt_mqs.end(), mqs.begin(), mqs.end());
  }

  return base::MannWhitneyEffectSize<u8>(ref_mqs, absl::MakeConstSpan(alt_mqs));
}

// ============================================================================
// ReadPosCohenD (RPCD FORMAT field)
//
// Coverage-normalized effect size comparing folded read positions of REF vs
// ALT reads. Folded position = min(p, 1−p) where p = variant_query_pos / read_len.
// Maps both 5' and 3' read edges to 0.0, centers to 0.5.
//
// True variants should be uniformly distributed across read positions.
// Artifacts from 3' quality degradation cluster at read edges (low folded
// position), producing a negative effect size for ALT.
//
// Returns 0.0 if either group is empty or all positions are identical.
// ============================================================================
auto VariantSupport::ReadPosCohenD() const -> f64 {
  absl::Span<const f64> ref_positions;
  if (REF_ALLELE_IDX < mAlleleData.size()) {
    ref_positions = absl::MakeConstSpan(mAlleleData[REF_ALLELE_IDX].folded_read_positions);
  }

  std::vector<f64> alt_positions;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    const auto& pos = mAlleleData[i].folded_read_positions;
    alt_positions.insert(alt_positions.end(), pos.begin(), pos.end());
  }

  return base::MannWhitneyEffectSize<f64>(ref_positions, absl::MakeConstSpan(alt_positions));
}

// ============================================================================
// BaseQualCohenD (BQCD FORMAT field)
//
// Coverage-normalized effect size comparing base qualities of REF vs ALT reads.
// Concatenates forward and reverse strand base qualities per allele.
// Detects 8-oxoguanine oxidation artifacts where the miscalled base has
// characteristically low Phred confidence.
//
// Returns 0.0 if either group is empty or all qualities are identical.
// ============================================================================
auto VariantSupport::BaseQualCohenD() const -> f64 {
  // Collect REF base qualities (fwd + rev concatenated)
  std::vector<u8> ref_bqs;
  if (REF_ALLELE_IDX < mAlleleData.size()) {
    const auto& ref = mAlleleData[REF_ALLELE_IDX];
    ref_bqs.reserve(ref.fwd_base_quals.size() + ref.rev_base_quals.size());
    ref_bqs.insert(ref_bqs.end(), ref.fwd_base_quals.begin(), ref.fwd_base_quals.end());
    ref_bqs.insert(ref_bqs.end(), ref.rev_base_quals.begin(), ref.rev_base_quals.end());
  }

  // Pool all ALT base qualities
  std::vector<u8> alt_bqs;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    const auto& alt = mAlleleData[i];
    alt_bqs.insert(alt_bqs.end(), alt.fwd_base_quals.begin(), alt.fwd_base_quals.end());
    alt_bqs.insert(alt_bqs.end(), alt.rev_base_quals.begin(), alt.rev_base_quals.end());
  }

  return base::MannWhitneyEffectSize<u8>(absl::MakeConstSpan(ref_bqs), absl::MakeConstSpan(alt_bqs));
}

// ============================================================================
// AlleleMismatchDelta (ASMD FORMAT field)
//
// mean(ALT NM) − mean(REF NM), where NM is the edit distance of each read
// against the REF haplotype. Both groups share the variant's own edit
// contribution, so ASMD cancels it and isolates excess noise.
//
// Positive ASMD → ALT reads have more mismatches (chimeric/paralogous signal)
// Near zero     → expected for true variants
// Returns 0.0 if either group is empty.
// ============================================================================
auto VariantSupport::AlleleMismatchDelta() const -> f64 {
  // REF group mean NM
  f64 ref_mean = 0.0;
  if (REF_ALLELE_IDX < mAlleleData.size()) {
    const auto& ref_nms = mAlleleData[REF_ALLELE_IDX].ref_nm_values;
    if (!ref_nms.empty()) {
      ref_mean = std::accumulate(ref_nms.begin(), ref_nms.end(), 0.0) / static_cast<f64>(ref_nms.size());
    }
  }

  // Pool all ALT NM values
  f64 alt_sum = 0.0;
  usize alt_count = 0;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    const auto& nms = mAlleleData[i].ref_nm_values;
    alt_sum += std::accumulate(nms.begin(), nms.end(), 0.0);
    alt_count += nms.size();
  }

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (alt_count == 0) return 0.0;
  const auto alt_mean = alt_sum / static_cast<f64>(alt_count);
  return alt_mean - ref_mean;
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
  // strand information is already captured in the SB (Phred-scaled Fisher strand bias) field.
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
