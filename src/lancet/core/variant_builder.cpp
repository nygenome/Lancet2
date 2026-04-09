#include "lancet/core/variant_builder.h"

#include "lancet/base/logging.h"
#include "lancet/base/repeat.h"
#include "lancet/base/sliding.h"
#include "lancet/base/types.h"
#include "lancet/caller/msa_builder.h"
#include "lancet/caller/variant_call.h"
#include "lancet/caller/variant_set.h"
#include "lancet/caller/variant_support.h"
#include "lancet/cbdg/read.h"
#include "lancet/core/sample_info.h"
#include "lancet/core/window.h"

#include "absl/hash/hash.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/core.h"
#include "spoa/alignment_engine.hpp"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <numeric>
#include <spdlog/fmt/bundled/format.h>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace lancet::core {

namespace {

/*
 * ============================================================================
 * SPOA MSA Parameter Rationale for Lancet2 Variant Extraction
 * ============================================================================
 * Values: Match: 0, Mismatch: -6, Gap1: -6,-2, Gap2: -26,-1
 *
 * **SIMD Overflow Note**: All classical parameters (typically +2/-4) have
 * been mathematically shifted downwards by 2. Setting the Match ceiling rigidly
 * to 0 is mandatory because SPOA employs 8-bit AVX2 SIMD vector lanes. A 1000bp
 * window accumulating +2 match scores would reach 2000, severely overflowing the
 * signed 8-bit integers (max 127). By anchoring the Match score to 0, all runtime
 * scores accumulate negatively, organically staying within the unrolled boundaries.
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
 *    mutation clusters. We use 0 / -6. The geometrical difference (6) keeps
 *    the MSA robustly intact while globally forcing alignments through complex variants.
 *
 * 3. Micro-Indel Sensitivity (Convex Model 1: -6, -2):
 *    asm5's -39 gap open penalty prevents small indels, forcing them to misalign
 *    as false-positive SNPs. Our -6 open / -2 extend penalty allows true small
 *    biological indels to open naturally while still applying enough friction
 *    to prevent 1bp sequencing errors (e.g., homopolymer stutters) from opening gaps.
 *
 * 4. Large Insertion/Deletion Continuity (Convex Model 2: -26, -1):
 *    asm5's -81 penalty for large gaps will soft-clip contigs right at an insertion/deletion
 *    breakpoint. Our parameters mathematically intersect at exactly 20bp (6 + 2L = 26 + 1L).
 *    For gaps > 20bp, the algorithm switches to Model 2 where the extension cost
 *    drops to -1. This "cheap extension" forces the DP matrix into mapping massive
 *    insertions/deletions as single, contiguous blocks in the MSA rather than
 *    dropping the alignment entirely.
 *
 *    https://curiouscoding.nl/posts/pairwise-alignment ->
 *    – Convex dual affine gap scoring -> min(g1+(i-1)*e1, g2+(i-1)*e2)
 */
constexpr i8 MSA_MATCH_SCORE = 0;
constexpr i8 MSA_MISMATCH_SCORE = -6;
constexpr i8 MSA_OPEN1_SCORE = -6;
constexpr i8 MSA_EXTEND1_SCORE = -2;
constexpr i8 MSA_OPEN2_SCORE = -26;
constexpr i8 MSA_EXTEND2_SCORE = -1;
constexpr u8 DNA_ALPHABET_SIZE = 4;
constexpr u32 PREALLOC_WINDOW_LENGTH_MULTIPLIER = 3;

}  // namespace

VariantBuilder::VariantBuilder(std::shared_ptr<Params const> params, u32 const window_length)
    : mDebruijnGraph(params->mGraphParams),
      mReadCollector(params->mRdCollParams),
      mParamsPtr(std::move(params)),
      mSpoaState{.mEngine = spoa::AlignmentEngine::Create(
                     spoa::AlignmentType::kNW, MSA_MATCH_SCORE, MSA_MISMATCH_SCORE, MSA_OPEN1_SCORE,
                     MSA_EXTEND1_SCORE, MSA_OPEN2_SCORE, MSA_EXTEND2_SCORE),
                 .mGraph = spoa::Graph()},
      mAnnotator(mParamsPtr->mGcFraction) {
  mSpoaState.mEngine->Prealloc(window_length * PREALLOC_WINDOW_LENGTH_MULTIPLIER,
                               DNA_ALPHABET_SIZE);
  mGenotyper.SetNumSamples(mParamsPtr->mRdCollParams.SamplesCount());
}

auto VariantBuilder::ProcessWindow(std::shared_ptr<Window const> const& window) -> WindowResults {
  auto const region = window->AsRegionPtr();
  auto const reg_str = region->ToSamtoolsRegion();
  static thread_local auto const THREAD_ID =
      absl::Hash<std::thread::id>()(std::this_thread::get_id());
  LOG_DEBUG("Processing window {} in thread {:#x}", reg_str, THREAD_ID)

  if (static_cast<usize>(std::ranges::count(window->SeqView(), 'N')) == window->Length()) {
    LOG_DEBUG("Skipping window {} since it has only N bases in reference", reg_str)
    mCurrentCode = StatusCode::SKIPPED_NONLY_REF_BASES;
    return {};
  }

  if (HasExactRepeat(SlidingView(window->SeqView(), mParamsPtr->mGraphParams.mMaxKmerLen))) {
    LOG_DEBUG("Skipping window {} since reference has repeat {}-mers", reg_str,
              mParamsPtr->mGraphParams.mMaxKmerLen)
    mCurrentCode = StatusCode::SKIPPED_REF_REPEAT_SEEN;
    return {};
  }

  auto const& rc_params = mParamsPtr->mRdCollParams;
  if (!mParamsPtr->mSkipActiveRegion && !ReadCollector::IsActiveRegion(rc_params, *region)) {
    LOG_DEBUG("Skipping window {} since it has no evidence of mutation in any sample", reg_str)
    mCurrentCode = StatusCode::SKIPPED_INACTIVE_REGION;
    return {};
  }

  LOG_DEBUG("Collecting all available sample reads for window {}", reg_str)
  auto const rc_result = mReadCollector.CollectRegionResult(*region);
  absl::Span<cbdg::Read const> const reads = absl::MakeConstSpan(rc_result.mSampleReads);
  absl::Span<SampleInfo const> const samples = absl::MakeConstSpan(rc_result.mSampleList);

  auto const total_cov = SampleInfo::CombinedSampledCov(samples, window->Length());
  if (total_cov < static_cast<f64>(mParamsPtr->mGraphParams.mMinAnchorCov)) {
    LOG_DEBUG("Skipping window {} since it has only {:.2f}x total sample coverage", reg_str,
              total_cov)
    mCurrentCode = StatusCode::SKIPPED_INACTIVE_REGION;
    return {};
  }

  LOG_DEBUG("Building graph for {} with {} sample reads and {:.2f}x total coverage", reg_str,
            reads.size(), total_cov)
  // First haplotype from each component will always be the reference haplotype sequence for the
  // graph
  auto const dbg_rslt = mDebruijnGraph.BuildComponentHaplotypes(window->AsRegionPtr(), reads);
  auto const& component_haplotypes = dbg_rslt.mGraphHaplotypes;

  static auto const SUMMER = [](u64 const sum, auto const& comp_haps) -> u64 {
    return sum + comp_haps.size() - 1;
  };
  auto const num_asm_haps =
      std::accumulate(component_haplotypes.cbegin(), component_haplotypes.cend(), 0, SUMMER);
  if (num_asm_haps == 0) {
    LOG_DEBUG("Could not assemble any haplotypes for window {} with k={}", reg_str,
              mDebruijnGraph.CurrentK())
    mCurrentCode = StatusCode::SKIPPED_NOASM_HAPLOTYPE;
    return {};
  }

  WindowResults variants;
  for (usize idx = 0; idx < component_haplotypes.size(); ++idx) {
    auto const nhaps = component_haplotypes[idx].size();
    auto const anchor_start = window->StartPos1() + dbg_rslt.mAnchorStartIdxs[idx];
    std::vector<cbdg::Graph::Path> const& comp_paths = component_haplotypes[idx];

    std::vector<std::string> comp_haps;
    comp_haps.reserve(comp_paths.size());
    for (auto const& path : comp_paths) {
      comp_haps.emplace_back(path.Sequence());
    }

    LOG_DEBUG("Building MSA for graph component {} from window {} with {} haplotypes", idx, reg_str,
              nhaps)

    absl::Span<std::string const> const ref_and_alt_haps = absl::MakeConstSpan(comp_haps);
    mSpoaState.UpdateSpoaState(ref_and_alt_haps);
    // only serialize if we have a path to serialize to
    mSpoaState.SerializeGraph(MakeGfaPath(*window, idx));

    caller::VariantSet vset;
    vset.ExtractVariantsFromGraph(mSpoaState.mGraph, *window, anchor_start);

    // Annotate complexity features if enabled — gated on CLI flags
    if (mParamsPtr->mEnableSequenceComplexity) {
      mAnnotator.AnnotateSequenceComplexity(vset, absl::MakeConstSpan(comp_haps));
    }
    if (mParamsPtr->mEnableGraphComplexity && idx < dbg_rslt.mComponentMetrics.size()) {
      VariantAnnotator::AnnotateGraphComplexity(vset, dbg_rslt.mComponentMetrics[idx]);
    }

    if (vset.IsEmpty()) {
      LOG_DEBUG("No variants found in graph component {} for window {} with {} haplotypes", idx,
                reg_str, nhaps)
      continue;
    }

    LOG_DEBUG("Found variant(s) in graph component {} for window {} with {} haplotypes", idx,
              reg_str, nhaps)
    auto genotyped = mGenotyper.Genotype(ref_and_alt_haps, reads, vset);

    // Drop variants where the graph assembled a haplotype but no reads actually support it.
    static auto const HAS_ALT_SUPPORT = [](caller::SupportArray const& evidence) -> bool {
      return std::ranges::any_of(evidence, [](auto const& item) -> bool {
        return item.mData && item.mData->TotalAltCov() > 0;
      });
    };

    std::ranges::for_each(vset, [&](auto const& var) -> void {
      caller::VariantCall::SupportsByVariant var_supports;
      if (auto iter = genotyped.find(&var);
          iter != genotyped.end() && HAS_ALT_SUPPORT(iter->second)) {
        var_supports.emplace(&var, std::move(iter->second));
      }

      if (var_supports.empty()) {
        return;
      }

      caller::VariantCall::FeatureFlags const feat_flags{
          .mEnableGraphComplexity = mParamsPtr->mEnableGraphComplexity,
          .mEnableSequenceComplexity = mParamsPtr->mEnableSequenceComplexity,
      };

      variants.emplace_back(std::make_unique<caller::VariantCall>(&var, std::move(var_supports),
                                                                  samples, feat_flags, total_cov));
    });
  }

  if (variants.empty()) {
    LOG_DEBUG("No variants found for window {} from {} assembled graph paths", reg_str,
              num_asm_haps)
    mCurrentCode = StatusCode::MISSING_NO_MSA_VARIANTS;
    return {};
  }

  mCurrentCode = StatusCode::FOUND_GENOTYPED_VARIANT;
  LOG_DEBUG("Genotyped {} variant(s) for window {} by re-aligning sample reads", variants.size(),
            reg_str)
  return variants;
}

auto VariantBuilder::MakeGfaPath(Window const& win, usize const comp_id) const
    -> std::filesystem::path {
  if (mParamsPtr->mOutGraphsDir.empty()) return {};

  auto const fname = fmt::format("msa__{}_{}_{}__c{}.gfa", win.ChromName(), win.StartPos1(),
                                 win.EndPos1(), comp_id);
  auto out_path = mParamsPtr->mOutGraphsDir / "poa_graph" / fname;
  std::filesystem::create_directories(mParamsPtr->mOutGraphsDir / "poa_graph");
  return out_path;
}

auto ToString(VariantBuilder::StatusCode const status_code) -> std::string {
  using enum VariantBuilder::StatusCode;

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
