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
  const auto haplotype_paths = mDebruijnGraph.MakeHaplotypes(window->AsRegionPtr(), reads);
  if (haplotype_paths.empty()) {
    LOG_DEBUG("Could not assemble any haplotypes for window {} with k={}", reg_str, mDebruijnGraph.CurrentK())
    mCurrentCode = StatusCode::SKIPPED_NOASM_HAPLOTYPE;
    return {};
  }

  const auto nhaps = haplotype_paths.size();
  LOG_DEBUG("Building POA based MSA for window {} with ref anchor and {} haplotypes", reg_str, nhaps)
  const absl::Span<const std::string> alt_seqs = absl::MakeConstSpan(haplotype_paths);
  const caller::MsaBuilder msa_builder(window->AsRegionPtr(), alt_seqs, MakeGfaPath(*window));
  const caller::VariantSet vset(msa_builder, *window);
  if (vset.IsEmpty()) {
    LOG_DEBUG("No variants found for window {} from MSA of ref anchor and {} haplotypes", reg_str, nhaps)
    mCurrentCode = StatusCode::MISSING_NO_MSA_VARIANTS;
    return {};
  }

  const auto var_cnt = vset.Count();
  LOG_DEBUG("Found {} variant(s) for window {} from MSA of ref anchor and {} haplotypes", var_cnt, reg_str, nhaps)

  const auto klen = mDebruijnGraph.CurrentK();
  const auto &vprms = mParamsPtr->mVariantParams;

  WindowResults variants;
  variants.reserve(var_cnt);

  for (auto &&[var_ptr, evidence] : mGenotyper.Genotype(*region, alt_seqs, reads, vset)) {
    variants.emplace_back(std::make_unique<caller::VariantCall>(var_ptr, std::move(evidence), samples, vprms, klen));
  }

  mCurrentCode = StatusCode::FOUND_GENOTYPED_VARIANT;
  LOG_DEBUG("Genotyped {} variant(s) for window {} by re-aligning sample reads", var_cnt, reg_str)
  return variants;
}

auto VariantBuilder::MakeGfaPath(const Window &window) const -> std::filesystem::path {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mParamsPtr->mOutGraphsDir.empty()) return {};

  const auto fname = fmt::format("msa__{}_{}_{}.gfa", window.ChromName(), window.StartPos1(), window.EndPos1());
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
