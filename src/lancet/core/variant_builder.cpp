#include "lancet/core/variant_builder.h"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <numeric>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "absl/hash/hash.h"
#include "absl/types/span.h"
#include "lancet/base/logging.h"
#include "lancet/base/repeat.h"
#include "lancet/base/sliding.h"
#include "lancet/base/types.h"
#include "lancet/caller/msa_builder.h"
#include "lancet/caller/variant_call.h"
#include "lancet/caller/variant_set.h"
#include "lancet/cbdg/read.h"
#include "lancet/core/sample_info.h"
#include "lancet/core/window.h"
#include "spdlog/fmt/bundled/core.h"
#include "spoa/alignment_engine.hpp"

namespace lancet::core {

namespace {

  /*
 * ============================================================================
 * SPOA MSA Parameter Rationale for Lancet2 Variant Extraction
 * ============================================================================
 * Values: Match: 2, Mismatch: -4, Gap1: -4,-2, Gap2: -24,-1
 *
 * Unlike minimap2's `asm5` preset (which aggressively splits contigs at major
 * divergences for whole-genome synteny filtering), these parameters are tuned
 * to force end-to-end global alignment within a specific micro-assembly window
 * to capture dense somatic mutations and large Insertions/Deletions.
 *
 * 1. Why Convex (Dual-Affine) vs. Affine or Linear Scoring:
 *    - Linear Scoring applies a flat penalty per gap base, which is biologically
 *      inaccurate (one 50bp deletion is one biological event, not fifty 1bp
 *      independent events).
 *    - Single Affine Scoring forces a compromise: tune for small variants (strict
 *      extension) and you penalize/clip large insertions/deletions; tune for large
 *      insertions/deletions (loose extension) and sequencer noise creates messy,
 *      spurious small gaps.
 *    - Convex (Dual-Affine) Scoring solves this by taking the minimum of two
 *      intersecting models. It is strict for short gaps to suppress sequencer
 *      noise, but switches to an incredibly cheap extension penalty for large
 *      biological gaps.
 *
 * 2. Mismatch Tolerance (Multi-Nucleotide Variants / MNVs):
 *    asm5 uses a +1 match / -19 mismatch, which shatters alignments at dense
 *    mutation clusters. We use +2 / -4, meaning only 2 matching bases are needed
 *    to offset a SNP. This keeps the MSA globally intact through complex variants.
 *
 * 3. Micro-Indel Sensitivity (Convex Model 1: -4, -2):
 *    asm5's -39 gap open penalty prevents small indels, forcing them to misalign
 *    as false-positive SNPs. Our -4 open / -2 extend penalty allows true small
 *    biological indels to open naturally while still applying enough friction
 *    to prevent 1bp sequencing errors (e.g., homopolymer stutters) from opening gaps.
 *
 * 4. Large Insertion/Deletion Continuity (Convex Model 2: -24, -1):
 *    asm5's -81 penalty for large gaps will soft-clip contigs right at an insertion/deletion
 *    breakpoint. Our parameters mathematically intersect at 20bp (4 + 2L = 24 + 1L).
 *    For gaps > 20bp, the algorithm switches to Model 2 where the extension cost
 *    drops to -1. This "cheap extension" forces the DP matrix into mapping massive
 *    insertions/deletions as single, contiguous blocks in the MSA rather than
 *    dropping the alignment.
 *
 *    https://curiouscoding.nl/posts/pairwise-alignment ->
 *    – Convex dual affine gap scoring -> min(g1+(i-1)*e1, g2+(i-1)*e2)
 */
constexpr i8 MSA_MATCH_SCORE = 2;
constexpr i8 MSA_MISMATCH_SCORE = -4;
constexpr i8 MSA_OPEN1_SCORE = -4;
constexpr i8 MSA_EXTEND1_SCORE = -2;
constexpr i8 MSA_OPEN2_SCORE = -24;
constexpr i8 MSA_EXTEND2_SCORE = -1;
constexpr u8 DNA_ALPHABET_SIZE = 4;
constexpr u32 PREALLOC_WINDOW_LENGTH_MULTIPLIER = 3;

}  // namespace

VariantBuilder::VariantBuilder(std::shared_ptr<const Params> params, const u32 window_length)
    : mDebruijnGraph(params->mGraphParams), mReadCollector(params->mRdCollParams), mParamsPtr(std::move(params)),
      mAlnEngine(spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, MSA_MATCH_SCORE, MSA_MISMATCH_SCORE,
        MSA_OPEN1_SCORE, MSA_EXTEND1_SCORE, MSA_OPEN2_SCORE, MSA_EXTEND2_SCORE)) {
  mAlnEngine->Prealloc(window_length * PREALLOC_WINDOW_LENGTH_MULTIPLIER, DNA_ALPHABET_SIZE);
  mGenotyper.SetNumSamples(mParamsPtr->mRdCollParams.SamplesCount());
}

auto VariantBuilder::ProcessWindow(const std::shared_ptr<const Window> &window) -> WindowResults {
  const auto region = window->AsRegionPtr();
  const auto reg_str = region->ToSamtoolsRegion();
  static thread_local const auto tid = absl::Hash<std::thread::id>()(std::this_thread::get_id());
  LOG_DEBUG("Processing window {} in thread {:#x}", reg_str, tid)

  if (static_cast<usize>(std::ranges::count(window->SeqView(), 'N')) == window->Length()) {
    LOG_DEBUG("Skipping window {} since it has only N bases in reference", reg_str)
    mCurrentCode = StatusCode::SKIPPED_NONLY_REF_BASES;
    return {};
  }

  if (HasExactRepeat(SlidingView(window->SeqView(), mParamsPtr->mGraphParams.mMaxKmerLen))) {
    LOG_DEBUG("Skipping window {} since reference has repeat {}-mers", reg_str, mParamsPtr->mGraphParams.mMaxKmerLen)
    mCurrentCode = StatusCode::SKIPPED_REF_REPEAT_SEEN;
    return {};
  }

  const auto &rc_params = mParamsPtr->mRdCollParams;
  if (!mParamsPtr->mSkipActiveRegion && !ReadCollector::IsActiveRegion(rc_params, *region)) {
    LOG_DEBUG("Skipping window {} since it has no evidence of mutation in any sample", reg_str)
    mCurrentCode = StatusCode::SKIPPED_INACTIVE_REGION;
    return {};
  }

  LOG_DEBUG("Collecting all available sample reads for window {}", reg_str)
  const auto rc_result = mReadCollector.CollectRegionResult(*region);
  const absl::Span<const cbdg::Read> reads = absl::MakeConstSpan(rc_result.mSampleReads);
  const absl::Span<const SampleInfo> samples = absl::MakeConstSpan(rc_result.mSampleList);

  const auto total_cov = SampleInfo::CombinedSampledCov(samples, window->Length());
  if (total_cov < static_cast<f64>(mParamsPtr->mGraphParams.mMinAnchorCov)) {
    LOG_DEBUG("Skipping window {} since it has only {:.2f}x total sample coverage", reg_str, total_cov)
    mCurrentCode = StatusCode::SKIPPED_INACTIVE_REGION;
    return {};
  }

  LOG_DEBUG("Building graph for {} with {} sample reads and {:.2f}x total coverage", reg_str, reads.size(), total_cov)
  // First haplotype from each component will always be the reference haplotype sequence for the graph
  const auto dbg_rslt = mDebruijnGraph.BuildComponentHaplotypes(window->AsRegionPtr(), reads);
  const auto &component_haplotypes = dbg_rslt.mGraphHaplotypes;

  static const auto summer = [](const u64 sum, const auto &comp_haps) -> u64 { return sum + comp_haps.size() - 1; };
  const auto num_asm_haps = std::accumulate(component_haplotypes.cbegin(), component_haplotypes.cend(), 0, summer);
  if (num_asm_haps == 0) {
    LOG_DEBUG("Could not assemble any haplotypes for window {} with k={}", reg_str, mDebruijnGraph.CurrentK())
    mCurrentCode = StatusCode::SKIPPED_NOASM_HAPLOTYPE;
    return {};
  }

  WindowResults variants;
  for (usize idx = 0; idx < component_haplotypes.size(); ++idx) {
    const auto nhaps = component_haplotypes[idx].size();
    const auto anchor_start = window->StartPos1() + dbg_rslt.mAnchorStartIdxs[idx];
    const std::vector<std::string> &comp_haps = component_haplotypes[idx];
    LOG_DEBUG("Building MSA for graph component {} from window {} with {} haplotypes", idx, reg_str, nhaps)

    const absl::Span<const std::string> ref_and_alt_haps = absl::MakeConstSpan(comp_haps);
    const caller::MsaBuilder msa_builder(ref_and_alt_haps, *mAlnEngine, MakeGfaPath(*window, idx));
    const caller::VariantSet vset(msa_builder, *window, anchor_start);

    ScoreVariantLCR(vset, absl::MakeConstSpan(comp_haps));

    if (vset.IsEmpty()) {
      LOG_DEBUG("No variants found in graph component {} for window {} with {} haplotypes", idx, reg_str, nhaps)
      continue;
    }

    LOG_DEBUG("Found variant(s) in graph component {} for window {} with {} haplotypes", idx, reg_str, nhaps)
    auto genotyped = mGenotyper.Genotype(ref_and_alt_haps, reads, vset);
    auto grouped = vset.GroupByLocus();

    // Drop variants where the graph assembled a haplotype but no reads actually support it.
    static const auto HasAltSupport = [](const caller::Genotyper::PerSampleEvidence& evidence) {
      return std::ranges::any_of(evidence, [](const auto& kv) { return kv.second && kv.second->TotalAltCov() > 0; });
    };

    for (const auto& [locus_key, locus_variants] : grouped) {
      caller::VariantCall::SupportsByVariant locus_supports;
      for (const auto* var_ptr : locus_variants) {
        auto it = genotyped.find(var_ptr);
        if (it == genotyped.end() || !HasAltSupport(it->second)) continue;
        locus_supports.emplace(var_ptr, std::move(it->second));
      }

      if (locus_supports.empty()) continue;
      variants.emplace_back(std::make_unique<caller::VariantCall>(
          absl::MakeConstSpan(locus_variants), std::move(locus_supports), samples,
          absl::MakeConstSpan(mParamsPtr->mAnnotationFeatures)));
    }
  }

  if (variants.empty()) {
    LOG_DEBUG("No variants found for window {} from {} assembled graph paths", reg_str, num_asm_haps)
    mCurrentCode = StatusCode::MISSING_NO_MSA_VARIANTS;
    return {};
  }

  mCurrentCode = StatusCode::FOUND_GENOTYPED_VARIANT;
  LOG_DEBUG("Genotyped {} variant(s) for window {} by re-aligning sample reads", variants.size(), reg_str)
  return variants;
}

// ============================================================================
// ScoreVariantLCR: multi-scale low-complexity scoring for ALT_LCR and REF_LCR
//
// For each variant in the set, computes longdust Q(x)/ℓ complexity scores
// at 5 flanking distances around the variant position in each haplotype.
//
//   ALT_LCR: centered on the variant position in each haplotype that carries
//            it (REF + ALTs). Takes the max score across all haplotypes.
//   REF_LCR: centered on the REF allele position in the reference haplotype
//            (haplotypes[0]) only.
//
//   Scale 0:   [---5bp----|VVVV|----5bp---]
//   Scale 1:   [--10bp----|VVVV|----10bp--]
//   Scale 2:   [--50bp----|VVVV|----50bp--]
//   Scale 3:   [-100bp----|VVVV|---100bp--]
//   Scale 4:   [entire haplotype sequence ]
//
// Gated on the --annotation-features parameter. If neither ALT_LCR nor
// REF_LCR is requested, this method is a no-op.
// ============================================================================
void VariantBuilder::ScoreVariantLCR(const caller::VariantSet& vset,
                                     absl::Span<const std::string> haplotypes) const {
  const auto& ann_features = mParamsPtr->mAnnotationFeatures;
  static const auto has_feature = [](absl::Span<const std::string> feats, std::string_view name) -> bool {
    return std::ranges::find(feats, name) != feats.end();
  };

  const bool compute_alt_lcr = has_feature(ann_features, "ALT_LCR");
  const bool compute_ref_lcr = has_feature(ann_features, "REF_LCR");
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (!compute_alt_lcr && !compute_ref_lcr) return;

  for (auto& var : vset) {
    auto alt_scores = base::DEFAULT_LCR_SCORES;
    auto ref_scores = base::DEFAULT_LCR_SCORES;

    // ── ALT_LCR: score across all haplotypes carrying this variant ──
    if (compute_alt_lcr) {
      const auto alt_var_len = static_cast<i64>(std::max(var.mRefAllele.size(), var.mAltAllele.size()));

      for (const auto& [hap_idx, hap_pos] : var.mHapStart0Idxs) {
        if (hap_idx >= haplotypes.size()) continue;
        const auto& hap = haplotypes[hap_idx];
        const auto hap_len = static_cast<i64>(hap.size());
        const auto pos = static_cast<i64>(hap_pos);

        for (usize scale = 0; scale < base::LCR_FLANKS.size(); ++scale) {
          const auto flank = base::LCR_FLANKS[scale];
          const auto start = std::max<i64>(0, pos - flank);
          const auto end = std::min<i64>(hap_len, pos + alt_var_len + flank);
          const auto subseq = std::string_view(hap).substr(start, end - start);
          alt_scores[scale] = std::max(alt_scores[scale], mFlankScorer.Score(subseq));
        }

        alt_scores[base::NUM_LCR_SCALES - 1] =
            std::max(alt_scores[base::NUM_LCR_SCALES - 1], mHaplotypeScorer.Score(hap));
      }
    }

    // ── REF_LCR: score only on reference haplotype (haplotypes[0]) ──
    if (compute_ref_lcr) {
      static constexpr usize REF_HAP_IDX = 0;
      const auto ref_pos_it = var.mHapStart0Idxs.find(REF_HAP_IDX);
      const auto ref_var_len = static_cast<i64>(var.mRefAllele.size());

      if (ref_pos_it != var.mHapStart0Idxs.end() && REF_HAP_IDX < haplotypes.size()) {
        const auto& ref_hap = haplotypes[REF_HAP_IDX];
        const auto ref_hap_len = static_cast<i64>(ref_hap.size());
        const auto ref_pos = static_cast<i64>(ref_pos_it->second);

        for (usize scale = 0; scale < base::LCR_FLANKS.size(); ++scale) {
          const auto flank = base::LCR_FLANKS[scale];
          const auto start = std::max<i64>(0, ref_pos - flank);
          const auto end = std::min<i64>(ref_hap_len, ref_pos + ref_var_len + flank);
          const auto subseq = std::string_view(ref_hap).substr(start, end - start);
          ref_scores[scale] = std::max(ref_scores[scale], mFlankScorer.Score(subseq));
        }

        ref_scores[base::NUM_LCR_SCALES - 1] = mHaplotypeScorer.Score(ref_hap);
      }
    }

    // mAltLcrScores/mRefLcrScores are not part of RawVariant ordering/hash, so const_cast is safe
    auto& mutable_var = const_cast<caller::RawVariant&>(var);
    mutable_var.mAltLcrScores = alt_scores;
    mutable_var.mRefLcrScores = ref_scores;
  }
}

auto VariantBuilder::MakeGfaPath(const Window &win, const usize comp_id) const -> std::filesystem::path {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mParamsPtr->mOutGraphsDir.empty()) return {};

  const auto fname = fmt::format("msa__{}_{}_{}__c{}.gfa", win.ChromName(), win.StartPos1(), win.EndPos1(), comp_id);
  auto out_path = mParamsPtr->mOutGraphsDir / "poa_graph" / fname;
  std::filesystem::create_directories(mParamsPtr->mOutGraphsDir / "poa_graph");
  return out_path;
}

auto ToString(const VariantBuilder::StatusCode status_code) -> std::string {
  using VariantBuilder::StatusCode::FOUND_GENOTYPED_VARIANT;
  using VariantBuilder::StatusCode::MISSING_NO_MSA_VARIANTS;
  using VariantBuilder::StatusCode::SKIPPED_INACTIVE_REGION;
  using VariantBuilder::StatusCode::SKIPPED_NOASM_HAPLOTYPE;
  using VariantBuilder::StatusCode::SKIPPED_NONLY_REF_BASES;
  using VariantBuilder::StatusCode::SKIPPED_REF_REPEAT_SEEN;

  switch (status_code) {
    case SKIPPED_NONLY_REF_BASES:
      return "SKIPPED_NONLY_REF_BASES";
    case SKIPPED_REF_REPEAT_SEEN:
      return "SKIPPED_REF_REPEAT_SEEN";
    case SKIPPED_INACTIVE_REGION:
      return "SKIPPED_INACTIVE_REGION";
    case SKIPPED_NOASM_HAPLOTYPE:
      return "SKIPPED_NOASM_HAPLOTYPE";
    case MISSING_NO_MSA_VARIANTS:
      return "MISSING_NO_MSA_VARIANTS";
    case FOUND_GENOTYPED_VARIANT:
      return "FOUND_GENOTYPED_VARIANT";
    default:
      break;
  }

  return "UNKNOWN";
}

}  // namespace lancet::core
