#include "lancet/caller/variant_support.h"

#include "lancet/base/mann_whitney.h"
#include "lancet/base/types.h"
#include "lancet/caller/genotype_likelihood.h"
#include "lancet/caller/posterior_base_qual.h"

#include "absl/types/span.h"

#include <algorithm>
#include <memory>
#include <numeric>
#include <string_view>
#include <utility>
#include <vector>

#include <cmath>

namespace lancet::caller {

// ============================================================================
// AddEvidence
// ============================================================================
void VariantSupport::AddEvidence(ReadEvidence const& evidence) {
  EnsureAlleleSlot(evidence.mAllele);
  auto& data = mAlleleData[evidence.mAllele];

  // Deduplicate by read name hash: only first-seen strand counts.
  auto const [iter, inserted] = data.mNameHashes.try_emplace(evidence.mRnameHash, evidence.mStrand);
  if (!inserted) return;

  // Each read contributes exactly one representative base quality value.
  // For indels, this is already collapsed to min(PBQ) across the variant region
  // by the genotyper (see ReadAlleleAssignment::mBaseQualAtVar comment).
  if (evidence.mStrand == Strand::FWD) {
    data.mFwdBaseQuals.push_back(evidence.mBaseQual);
  } else {
    data.mRevBaseQuals.push_back(evidence.mBaseQual);
  }

  data.mMapQuals.push_back(evidence.mMapQual);
  data.mAlnScores.push_back(evidence.mAlnScore);

  // Track soft-clip status from original alignment (for SCA FORMAT tag)
  if (evidence.mIsSoftClipped) ++data.mSoftClipCount;

  // Track insert sizes from properly-paired reads (for FLD FORMAT tag).
  // Signed insert size preserved: negative = ALT shorter (cfDNA biology).
  if (evidence.mIsProperPair && evidence.mInsertSize != 0) {
    data.mProperPairIsizes.push_back(static_cast<f64>(evidence.mInsertSize));
  }

  // Track folded read position (for RPCD FORMAT tag)
  data.mFoldedReadPositions.push_back(evidence.mFoldedReadPos);

  // Track edit distance to REF haplotype (for ASMD FORMAT tag)
  data.mRefNmValues.push_back(static_cast<f64>(evidence.mRefNm));

  // Track fragment start position (for FSSE: PCR optical duplicate detection)
  data.mAlignmentStarts.push_back(evidence.mAlignmentStart);

  // Track edit distance against assigned haplotype (for AHDD)
  data.mOwnHapNmValues.push_back(static_cast<f64>(evidence.mOwnHapNm));

  // Track assigned SPOA haplotype ID (for HSE: path co-segregation)
  data.mHaplotypeIds.push_back(evidence.mAssignedHaplotypeId);
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
void VariantSupport::MergeAlleleFrom(VariantSupport const& src, AlleleIndex const src_allele,
                                     AlleleIndex const dst_allele) {
  if (src_allele >= src.mAlleleData.size()) return;

  EnsureAlleleSlot(dst_allele);
  auto const& src_data = src.mAlleleData[src_allele];
  auto& dst_data = mAlleleData[dst_allele];

  for (auto const& [rname_hash, strand] : src_data.mNameHashes) {
    auto const [iter, inserted] = dst_data.mNameHashes.try_emplace(rname_hash, strand);
    // already seen this read in dst
    if (!inserted) continue;
  }

  // Append PBQ, RMQ, and alignment scores (reads already deduped above by mNameHashes,
  // but the qual vectors are append-only during AddEvidence, so we append here too).
  dst_data.mFwdBaseQuals.insert(dst_data.mFwdBaseQuals.end(), src_data.mFwdBaseQuals.begin(),
                                src_data.mFwdBaseQuals.end());

  dst_data.mRevBaseQuals.insert(dst_data.mRevBaseQuals.end(), src_data.mRevBaseQuals.begin(),
                                src_data.mRevBaseQuals.end());

  dst_data.mMapQuals.insert(dst_data.mMapQuals.end(), src_data.mMapQuals.begin(),
                            src_data.mMapQuals.end());

  dst_data.mAlnScores.insert(dst_data.mAlnScores.end(), src_data.mAlnScores.begin(),
                             src_data.mAlnScores.end());

  dst_data.mSoftClipCount += src_data.mSoftClipCount;

  dst_data.mProperPairIsizes.insert(dst_data.mProperPairIsizes.end(),
                                    src_data.mProperPairIsizes.begin(),
                                    src_data.mProperPairIsizes.end());

  dst_data.mFoldedReadPositions.insert(dst_data.mFoldedReadPositions.end(),
                                       src_data.mFoldedReadPositions.begin(),
                                       src_data.mFoldedReadPositions.end());

  dst_data.mRefNmValues.insert(dst_data.mRefNmValues.end(), src_data.mRefNmValues.begin(),
                               src_data.mRefNmValues.end());

  // Merge new artifact metric vectors
  dst_data.mAlignmentStarts.insert(dst_data.mAlignmentStarts.end(),
                                   src_data.mAlignmentStarts.begin(),
                                   src_data.mAlignmentStarts.end());

  dst_data.mOwnHapNmValues.insert(dst_data.mOwnHapNmValues.end(), src_data.mOwnHapNmValues.begin(),
                                  src_data.mOwnHapNmValues.end());

  dst_data.mHaplotypeIds.insert(dst_data.mHaplotypeIds.end(), src_data.mHaplotypeIds.begin(),
                                src_data.mHaplotypeIds.end());
}

// ============================================================================
// Per-Allele Accessors
// ============================================================================
auto VariantSupport::FwdCount(AlleleIndex const idx) const -> usize {
  if (idx >= mAlleleData.size()) return 0;
  return mAlleleData[idx].mFwdBaseQuals.size();
}

auto VariantSupport::RevCount(AlleleIndex const idx) const -> usize {
  if (idx >= mAlleleData.size()) return 0;
  return mAlleleData[idx].mRevBaseQuals.size();
}

auto VariantSupport::TotalAlleleCov(AlleleIndex const idx) const -> usize {
  return FwdCount(idx) + RevCount(idx);
}

auto VariantSupport::TotalSampleCov() const noexcept -> usize {
  return std::transform_reduce(
      mAlleleData.cbegin(), mAlleleData.cend(), usize{0}, std::plus<>{},
      [](auto const& entry) { return entry.mFwdBaseQuals.size() + entry.mRevBaseQuals.size(); });
}

auto VariantSupport::TotalAltCov() const -> usize {
  usize total = 0;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    total += TotalAlleleCov(static_cast<AlleleIndex>(i));
  }
  return total;
}

// ============================================================================
// RawPosteriorBaseQual — delegates to posterior_base_qual.h/.cpp
// See posterior_base_qual.h for Edgar & Flyvbjerg (2014) derivation.
// ============================================================================
auto VariantSupport::RawPosteriorBaseQual(AlleleIndex const idx) const -> f64 {
  if (idx >= mAlleleData.size()) return 0.0;
  auto const& data = mAlleleData[idx];
  return ComputeRawPosteriorBaseQual(absl::MakeConstSpan(data.mFwdBaseQuals),
                                     absl::MakeConstSpan(data.mRevBaseQuals));
}

// ============================================================================
// RmsMappingQual (RMQ FORMAT field)
//
// Root-mean-square of the mapping qualities of all reads assigned to this
// allele. This follows the standard INFO/RMQ convention (samtools, GATK).
//
//   RMQ = sqrt( Σ(mapq_i²) / N )
// ============================================================================
auto VariantSupport::RmsMappingQual(AlleleIndex const idx) const -> f64 {
  if (idx >= mAlleleData.size() || mAlleleData[idx].mMapQuals.empty()) return 0.0;

  auto const& mqs = mAlleleData[idx].mMapQuals;
  auto const square_summer = [](f64 const acc, u8 const mapq) -> f64 {
    return acc + (static_cast<f64>(mapq) * static_cast<f64>(mapq));
  };

  f64 const sum_sq = std::accumulate(mqs.cbegin(), mqs.cend(), 0.0, square_summer);
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
  auto const ref_fwd = static_cast<int>(FwdCount(REF_ALLELE_IDX));
  auto const ref_rev = static_cast<int>(RevCount(REF_ALLELE_IDX));

  // ALT strand counts (summed across all non-REF alleles)
  int alt_fwd = 0;
  int alt_rev = 0;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    alt_fwd += static_cast<int>(mAlleleData[i].mFwdBaseQuals.size());
    alt_rev += static_cast<int>(mAlleleData[i].mRevBaseQuals.size());
  }

  // Haldane correction: +1 to all cells handles zero-count edge cases
  // without sentinels. When all counts are zero, returns 0.0 (no bias).
  auto const rf1 = static_cast<f64>(ref_fwd + 1);
  auto const rr1 = static_cast<f64>(ref_rev + 1);
  auto const af1 = static_cast<f64>(alt_fwd + 1);
  auto const ar1 = static_cast<f64>(alt_rev + 1);
  return std::log((rf1 * ar1) / (rr1 * af1));
}

// ============================================================================
// MeanAlnScore
// ============================================================================
auto VariantSupport::MeanAlnScore(AlleleIndex const idx) const -> f64 {
  if (idx >= mAlleleData.size() || mAlleleData[idx].mAlnScores.empty()) return 0.0;

  auto const& scores = mAlleleData[idx].mAlnScores;
  f64 const sum = std::accumulate(scores.cbegin(), scores.cend(), 0.0);
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
    alt_sc += mAlleleData[i].mSoftClipCount;
    alt_total += TotalAlleleCov(static_cast<AlleleIndex>(i));
  }

  // REF soft-clip fraction
  auto const has_ref = REF_ALLELE_IDX < mAlleleData.size();
  auto const ref_sc = has_ref ? mAlleleData[REF_ALLELE_IDX].mSoftClipCount : usize{0};
  auto const ref_total = TotalAlleleCov(REF_ALLELE_IDX);

  f64 const alt_frac = alt_total > 0 ? static_cast<f64>(alt_sc) / static_cast<f64>(alt_total) : 0.0;
  f64 const ref_frac = ref_total > 0 ? static_cast<f64>(ref_sc) / static_cast<f64>(ref_total) : 0.0;
  return alt_frac - ref_frac;
}

// ============================================================================
// FragLengthDelta (FLD FORMAT field) — delegates to MeanAltMinusRef.
// See variant_support.h for the full derivation.
// ============================================================================
auto VariantSupport::FragLengthDelta() const -> std::optional<f64> {
  return MeanAltMinusRef<f64>(
      [](auto const& allele) -> auto const& { return allele.mProperPairIsizes; });
}

// ============================================================================
// MappingQualCohenD (MQCD FORMAT field) — delegates to RefVsAltEffectSize.
// See variant_support.h for the full derivation.
// ============================================================================
auto VariantSupport::MappingQualCohenD() const -> std::optional<f64> {
  return RefVsAltEffectSize<u8>([](auto const& allele) -> auto const& { return allele.mMapQuals; });
}

// ============================================================================
// ReadPosCohenD (RPCD FORMAT field) — delegates to RefVsAltEffectSize.
// See variant_support.h for the full derivation.
// ============================================================================
auto VariantSupport::ReadPosCohenD() const -> std::optional<f64> {
  return RefVsAltEffectSize<f64>(
      [](auto const& allele) -> auto const& { return allele.mFoldedReadPositions; });
}

// ============================================================================
// BaseQualCohenD (BQCD FORMAT field)
//
// Coverage-normalized effect size comparing base qualities of REF vs ALT reads.
// Concatenates forward and reverse strand base qualities per allele.
// Detects 8-oxoguanine oxidation artifacts where the miscalled base has
// characteristically low Phred confidence.
//
// Custom shape: fwd+rev concatenation per allele prevents using the simple
// RefVsAltEffectSize template (which expects a single field per allele).
//
// Returns std::nullopt if either group is empty (untestable).
// Returns 0.0 when test ran but found no bias (genuine zero).
// ============================================================================
auto VariantSupport::BaseQualCohenD() const -> std::optional<f64> {
  // Collect REF base qualities (fwd + rev concatenated)
  std::vector<u8> ref_bqs;
  if (REF_ALLELE_IDX < mAlleleData.size()) {
    auto const& ref = mAlleleData[REF_ALLELE_IDX];
    ref_bqs.reserve(ref.mFwdBaseQuals.size() + ref.mRevBaseQuals.size());
    ref_bqs.insert(ref_bqs.end(), ref.mFwdBaseQuals.begin(), ref.mFwdBaseQuals.end());
    ref_bqs.insert(ref_bqs.end(), ref.mRevBaseQuals.begin(), ref.mRevBaseQuals.end());
  }

  // Pool all ALT base qualities
  std::vector<u8> alt_bqs;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    auto const& alt = mAlleleData[i];
    alt_bqs.insert(alt_bqs.end(), alt.mFwdBaseQuals.begin(), alt.mFwdBaseQuals.end());
    alt_bqs.insert(alt_bqs.end(), alt.mRevBaseQuals.begin(), alt.mRevBaseQuals.end());
  }

  return base::MannWhitneyEffectSize<u8>(absl::MakeConstSpan(ref_bqs),
                                         absl::MakeConstSpan(alt_bqs));
}

// ============================================================================
// AlleleMismatchDelta (ASMD FORMAT field) — delegates to MeanAltMinusRef.
// Subtracts variant_length to isolate excess noise beyond the expected
// structural difference. See variant_support.h for the full derivation.
// ============================================================================
auto VariantSupport::AlleleMismatchDelta(usize const variant_length) const -> std::optional<f64> {
  return MeanAltMinusRef<f64>([](auto const& allele) -> auto const& { return allele.mRefNmValues; },
                              static_cast<f64>(variant_length));
}

// ============================================================================
// ComputeFSSE (FSSE FORMAT field) — delegates to AltPooledEntropy.
// 3 bp binning absorbs exonuclease fraying. Cap at 20 bins.
// See variant_support.h for the full derivation.
// ============================================================================
auto VariantSupport::ComputeFSSE() const -> std::optional<f64> {
  return AltPooledEntropy<i64, i64>(
      [](auto const& allele) -> auto const& { return allele.mAlignmentStarts; },
      [](i64 start) { return start / 3; },  // 3bp binning for exonuclease fraying
      20.0);                                // cap at 20 bins
}

// ============================================================================
// ComputeAHDD (AHDD FORMAT field) — delegates to MeanAltMinusRef.
// See variant_support.h for the full derivation.
// ============================================================================
auto VariantSupport::ComputeAHDD() const -> std::optional<f64> {
  return MeanAltMinusRef<f64>(
      [](auto const& allele) -> auto const& { return allele.mOwnHapNmValues; });
}

// ============================================================================
// ComputeHSE (HSE FORMAT field) — delegates to AltPooledEntropy.
// Normalized by log₂(total_haplotypes).
// See variant_support.h for the full derivation.
// ============================================================================
auto VariantSupport::ComputeHSE(usize const total_haplotypes) const -> std::optional<f64> {
  if (total_haplotypes < 2) return std::nullopt;
  return AltPooledEntropy<u32, u32>(
      [](auto const& allele) -> auto const& { return allele.mHaplotypeIds; },
      [](u32 hapid) { return hapid; }, static_cast<f64>(total_haplotypes));
}

// ============================================================================
// Genotype Likelihood Architecture
//
// The DM probability engine and CMLOD mixture model have been extracted to:
//   genotype_likelihood.h/.cpp
//
// ┌─────────────────────────────────────────────┐
// │ genotype_likelihood.h/.cpp                  │
// │   ComputeGenotypePLs(allele_counts)         │ ← DM model
// │   ComputeGenotypeQuality(pls)               │ ← GQ from PLs
// │   ComputeContinuousMixtureLods(quals, covs) │ ← CMLOD engine
// │   AlleleBaseQuals struct                    │
// └──────────────────────┬──────────────────────┘
//                        │ called by
// ┌──────────────────────▼──────────────────────┐
// │ variant_support.cpp (this file)             │
// │   ComputePLs()      → thin wrapper          │
// │   ComputeGQ()       → thin wrapper          │
// │   ComputeContinuousMixtureLods() → wrapper  │
// └─────────────────────────────────────────────┘
// ============================================================================

auto VariantSupport::ComputePLs(usize const num_alleles) const -> absl::InlinedVector<u32, 6> {
  auto const num_al = static_cast<int>(num_alleles);
  if (num_al == 0) return {};

  // Build allele count vector: count[i] = total reads assigned to allele i.
  // TotalAlleleCov(idx) returns 0 for alleles beyond mAlleleData.size(),
  // so alleles the sample has no reads for correctly get count=0.
  std::vector<int> allele_counts(num_al, 0);
  for (int allele_idx = 0; allele_idx < num_al; ++allele_idx) {
    auto const aidx = static_cast<AlleleIndex>(allele_idx);
    allele_counts[allele_idx] = static_cast<int>(TotalAlleleCov(aidx));
  }

  return ComputeGenotypePLs(absl::MakeConstSpan(allele_counts));
}

auto VariantSupport::ComputeGQ(absl::Span<u32 const> phred_likelihoods) -> u32 {
  return ComputeGenotypeQuality(phred_likelihoods);
}

auto VariantSupport::ComputeContinuousMixtureLods(usize const num_alleles) const
    -> std::vector<f64> {
  auto const num_al = static_cast<int>(num_alleles);
  // Braced init `{num_al, 0.0}` would invoke vector's initializer_list<double> ctor
  // (always size 2), silently losing the size-based fill semantics: num_al==0 must
  // produce an empty vector, num_al==1 must produce {0.0}. Required by the per-allele
  // CMLOD contract (REF at index 0 = 0.0; one score per allele).
  // NOLINTNEXTLINE(modernize-return-braced-init-list)
  if (num_al < 2) return std::vector<f64>(num_al, 0.0);

  // Build AlleleBaseQuals views into private PerAlleleData.
  // For alleles beyond mAlleleData.size(): empty qual spans + 0 coverage.
  // This produces CMLOD=0.0 for unobserved alleles (no evidence either way).
  std::vector<AlleleBaseQuals> allele_quals(num_al);
  std::vector<usize> allele_coverages(num_al);
  for (int index = 0; index < num_al; ++index) {
    auto const aidx = static_cast<usize>(index);
    if (aidx < mAlleleData.size()) {
      auto const& per_allele = mAlleleData[aidx];
      allele_quals[index] = AlleleBaseQuals{
          .mFwdBaseQuals = absl::MakeConstSpan(per_allele.mFwdBaseQuals),
          .mRevBaseQuals = absl::MakeConstSpan(per_allele.mRevBaseQuals),
      };
    }
    allele_coverages[index] = TotalAlleleCov(static_cast<AlleleIndex>(index));
  }

  return caller::ComputeContinuousMixtureLods(absl::MakeConstSpan(allele_quals),
                                              absl::MakeConstSpan(allele_coverages));
}

// ============================================================================
// EnsureAlleleSlot
// ============================================================================
void VariantSupport::EnsureAlleleSlot(AlleleIndex const idx) {
  if (idx >= mAlleleData.size()) mAlleleData.resize(idx + 1);
}

}  // namespace lancet::caller
