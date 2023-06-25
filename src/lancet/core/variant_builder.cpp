#include "lancet/core/variant_builder.h"

#include <thread>
#include <utility>

#include "absl/hash/hash.h"
#include "lancet/base/logging.h"
#include "lancet/base/repeat.h"
#include "lancet/base/sliding.h"
#include "lancet/caller/msa_builder.h"
#include "lancet/caller/variant_set.h"

namespace lancet::core {

VariantBuilder::VariantBuilder(std::shared_ptr<const Params> params)
    : mDebruijnGraph(params->mGraphParams), mReadCollector(params->mRdCollParams), mParamsPtr(std::move(params)) {
  mGenotyper.SetNumSamples(mParamsPtr->mRdCollParams.SamplesCount());
}

auto VariantBuilder::ProcessWindow(const std::shared_ptr<const Window> &window) -> WindowResults {
  const auto region = window->AsRegionPtr();
  const auto reg_str = region->ToSamtoolsRegion();
  static thread_local const auto tid = absl::Hash<std::thread::id>()(std::this_thread::get_id());
  LOG_DEBUG("Starting to process window {} in thread {:#x}", reg_str, tid)

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

  const auto rc_result = mReadCollector.CollectRegionResult(*region);
  const absl::Span<const cbdg::Read> reads = absl::MakeConstSpan(rc_result.mSampleReads);
  const absl::Span<const SampleInfo> samples = absl::MakeConstSpan(rc_result.mSampleList);

  const auto total_cov = SampleInfo::TotalMeanCov(samples, window->Length());
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
  const auto dbg_klen = mDebruijnGraph.CurrentK();
  const auto &vprms = mParamsPtr->mVariantParams;

  for (usize idx = 0; idx < component_haplotypes.size(); ++idx) {
    const auto nseqs = component_haplotypes[idx].size();
    const auto anchor_start = window->StartPos1() + dbg_rslt.mAnchorStartIdxs[idx];
    const std::vector<std::string> &comp_haps = component_haplotypes[idx];
    LOG_DEBUG("Building MSA for graph component {} from window {} with {} sequences", idx, reg_str, nseqs)

    const absl::Span<const std::string> ref_and_alt_haps = absl::MakeConstSpan(comp_haps);
    const caller::MsaBuilder msa_builder(ref_and_alt_haps, MakeGfaPath(*window, idx));
    const caller::VariantSet vset(msa_builder, *window, anchor_start);

    if (vset.IsEmpty()) {
      LOG_DEBUG("No variants found in graph component {} from window {} with {} sequences", idx, reg_str, nseqs)
      continue;
    }

    LOG_DEBUG("Found variant(s) in graph component {} from window {} with {} sequences", idx, reg_str, nseqs)
    for (auto &&[var, evidence] : mGenotyper.Genotype(ref_and_alt_haps, reads, vset)) {
      variants.emplace_back(std::make_unique<caller::VariantCall>(var, std::move(evidence), samples, vprms, dbg_klen));
    }
  }

  if (variants.empty()) {
    LOG_DEBUG("No variants found for window {} from {} assembled graph paths", reg_str, num_asm_haps)
    mCurrentCode = StatusCode::MISSING_NO_MSA_VARIANTS;
    return {};
  }

  mCurrentCode = StatusCode::FOUND_GENOTYPED_VARIANT;
  LOG_DEBUG("Genotyped variant(s) for window {} by re-aligning sample reads", reg_str)
  return variants;
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
