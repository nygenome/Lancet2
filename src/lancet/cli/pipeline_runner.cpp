#include "lancet/cli/pipeline_runner.h"

#include <algorithm>
#include <chrono>  // NOLINT(misc-include-cleaner)
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <memory>
#include <numeric>
#include <stop_token>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#ifdef LANCET_PROFILE_MODE
#include "gperftools/profiler.h"
#endif

#include "absl/container/btree_map.h"
#include "absl/container/fixed_array.h"
#include "absl/hash/hash.h"
#include "absl/strings/match.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "absl/time/clock.h"
#include "absl/time/time.h"
#include "absl/types/span.h"
#include "concurrentqueue.h"
#include "lancet/base/logging.h"
#include "lancet/base/timer.h"
#include "lancet/base/types.h"
#include "lancet/base/version.h"
#include "lancet/cli/cli_params.h"
#include "lancet/cli/eta_timer.h"
#include "lancet/core/async_worker.h"
#include "lancet/core/read_collector.h"
#include "lancet/core/variant_builder.h"
#include "lancet/core/variant_store.h"
#include "lancet/core/window.h"
#include "lancet/core/window_builder.h"
#include "lancet/hts/alignment.h"
#include "lancet/hts/bgzf_ostream.h"
#include "lancet/hts/extractor.h"
#include "lancet/hts/reference.h"
#include "spdlog/fmt/bundled/core.h"

using lancet::core::AsyncWorker;
using lancet::core::VariantBuilder;
using WindowStats = absl::btree_map<VariantBuilder::StatusCode, u64>;

namespace {

[[nodiscard]] inline auto InitWindowStats() -> absl::btree_map<VariantBuilder::StatusCode, u64> {
  using VariantBuilder::StatusCode::FOUND_GENOTYPED_VARIANT;
  using VariantBuilder::StatusCode::MISSING_NO_MSA_VARIANTS;
  using VariantBuilder::StatusCode::SKIPPED_INACTIVE_REGION;
  using VariantBuilder::StatusCode::SKIPPED_NOASM_HAPLOTYPE;
  using VariantBuilder::StatusCode::SKIPPED_NONLY_REF_BASES;
  using VariantBuilder::StatusCode::SKIPPED_REF_REPEAT_SEEN;
  using VariantBuilder::StatusCode::UNKNOWN;

  return WindowStats{{UNKNOWN, 0},
                     {SKIPPED_NONLY_REF_BASES, 0},
                     {SKIPPED_REF_REPEAT_SEEN, 0},
                     {SKIPPED_INACTIVE_REGION, 0},
                     {SKIPPED_NOASM_HAPLOTYPE, 0},
                     {MISSING_NO_MSA_VARIANTS, 0},
                     {FOUND_GENOTYPED_VARIANT, 0}};
}

void LogWindowStats(const WindowStats &stats) {
  using VariantBuilder::StatusCode::FOUND_GENOTYPED_VARIANT;
  using VariantBuilder::StatusCode::MISSING_NO_MSA_VARIANTS;
  using VariantBuilder::StatusCode::SKIPPED_INACTIVE_REGION;
  using VariantBuilder::StatusCode::SKIPPED_NOASM_HAPLOTYPE;
  using VariantBuilder::StatusCode::SKIPPED_NONLY_REF_BASES;
  using VariantBuilder::StatusCode::SKIPPED_REF_REPEAT_SEEN;

  using CodeCounts = std::pair<const VariantBuilder::StatusCode, u64>;
  static const auto summer = [](const u64 sum, const CodeCounts &item) -> u64 { return sum + item.second; };
  const auto nwindows = std::accumulate(stats.cbegin(), stats.cend(), 0, summer);

  std::ranges::for_each(stats, [&nwindows](const CodeCounts &item) {
    const auto [status_code, count] = item;
    const auto pct_count = (100.0 * static_cast<f64>(count)) / static_cast<f64>(nwindows);
    switch (status_code) {
      case SKIPPED_NONLY_REF_BASES:
        LOG_INFO("SKIPPED_NONLY_REF_BASES | {:>8.4f}% of total windows | {} windows", pct_count, count)
        break;
      case SKIPPED_REF_REPEAT_SEEN:
        LOG_INFO("SKIPPED_REF_REPEAT_SEEN | {:>8.4f}% of total windows | {} windows", pct_count, count)
        break;
      case SKIPPED_INACTIVE_REGION:
        LOG_INFO("SKIPPED_INACTIVE_REGION | {:>8.4f}% of total windows | {} windows", pct_count, count)
        break;
      case SKIPPED_NOASM_HAPLOTYPE:
        LOG_INFO("SKIPPED_NOASM_HAPLOTYPE | {:>8.4f}% of total windows | {} windows", pct_count, count)
        break;
      case MISSING_NO_MSA_VARIANTS:
        LOG_INFO("MISSING_NO_MSA_VARIANTS | {:>8.4f}% of total windows | {} windows", pct_count, count)
        break;
      case FOUND_GENOTYPED_VARIANT:
        LOG_INFO("FOUND_GENOTYPED_VARIANT | {:>8.4f}% of total windows | {} windows", pct_count, count)
        break;
      default:
        break;
    }
  });
}

}  // namespace

// NOLINTBEGIN(bugprone-easily-swappable-parameters,performance-unnecessary-value-param)
void PipelineWorker(std::stop_token stop_token, const moodycamel::ProducerToken *in_token,
                    AsyncWorker::InQueuePtr in_queue, AsyncWorker::OutQueuePtr out_queue,
                    AsyncWorker::VariantStorePtr vstore, AsyncWorker::BuilderParamsPtr params, u32 window_length) {
  // NOLINTEND(bugprone-easily-swappable-parameters,performance-unnecessary-value-param)
#ifdef LANCET_PROFILE_MODE
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (ProfilingIsEnabledForAllThreads() != 0) ProfilerRegisterThread();
#endif
  auto worker = std::make_unique<AsyncWorker>(std::move(in_queue), std::move(out_queue), std::move(vstore),
                                              std::move(params), window_length);
  worker->Process(std::move(stop_token), *in_token);
}

namespace lancet::cli {

PipelineRunner::PipelineRunner(std::shared_ptr<CliParams> params) : mParamsPtr(std::move(params)) {
#ifdef LANCET_PROFILE_MODE
  setenv("CPUPROFILE_PER_THREAD_TIMERS", "1", 1);
  setenv("CPUPROFILE_FREQUENCY", "10000", 1);
  const auto timestamp = absl::FormatTime("%Y%m%d%ET%H%M%S", absl::Now(), absl::LocalTimeZone());
  const auto fname = fmt::format("Lancet.cpu_profile.{}.bin", timestamp);
  ProfilerStart(fname.c_str());
#endif
}

// ---------------------------------------------------------------------------
// InitWindowBuilder: setup and sort regions for batch emission
// ---------------------------------------------------------------------------

auto PipelineRunner::InitWindowBuilder(const CliParams &params) -> core::WindowBuilder {
  core::WindowBuilder window_builder(params.mVariantBuilder.mRdCollParams.mRefPath, params.mWindowBuilder);
  window_builder.AddBatchRegions(absl::MakeConstSpan(params.mInRegions));
  window_builder.AddBatchRegions(params.mBedFile);

  if (window_builder.IsEmpty()) {
    LOG_WARN("No input regions provided to build windows. Using contigs in reference as input regions")
    window_builder.AddAllReferenceRegions();
  }

  // Sort input regions before batch emission to ensure deterministic genomic ordering
  window_builder.SortInputRegions();
  return window_builder;
}

// ---------------------------------------------------------------------------
// Run: primary pipeline execution with dynamic batch-fed window queue
// ---------------------------------------------------------------------------

void PipelineRunner::Run() {
  Timer timer;
  static thread_local const auto tid = std::this_thread::get_id();
  LOG_INFO("Using main thread {:#x} to synchronize variant calling pipeline", absl::Hash<std::thread::id>()(tid))

  ValidateAndPopulateParams();
  if (!mParamsPtr->mVariantBuilder.mOutGraphsDir.empty()) {
    mParamsPtr->mVariantBuilder.mGraphParams.mOutGraphsDir = mParamsPtr->mVariantBuilder.mOutGraphsDir;
    std::filesystem::remove_all(mParamsPtr->mVariantBuilder.mOutGraphsDir);
    std::filesystem::create_directories(mParamsPtr->mVariantBuilder.mOutGraphsDir);
  }

  mParamsPtr->mOutVcfGz = std::filesystem::absolute(mParamsPtr->mOutVcfGz);
  if (!std::filesystem::exists(mParamsPtr->mOutVcfGz.parent_path())) {
    std::filesystem::create_directories(mParamsPtr->mOutVcfGz.parent_path());
  }

  hts::BgzfOstream output_vcf;
  if (!output_vcf.Open(mParamsPtr->mOutVcfGz, hts::BgzfFormat::VCF)) {
    LOG_CRITICAL("Could not open output VCF file: {}", mParamsPtr->mOutVcfGz.string())
    std::exit(EXIT_FAILURE);
  }

  output_vcf << BuildVcfHeader(*mParamsPtr);

  // Initialize the window builder with sorted regions
  auto window_builder = InitWindowBuilder(*mParamsPtr);
  const auto num_total_windows = window_builder.ExpectedTargetWindows();
  LOG_INFO("Processing {} window(s) with {} VariantBuilder thread(s)", num_total_windows, mParamsPtr->mNumWorkerThreads)

  // Use the BATCH_SIZE from WindowBuilder for queue feeding granularity
  static constexpr auto BATCH_SIZE = core::WindowBuilder::BATCH_SIZE;

  // When total windows fit within a multiple of BATCH_SIZE, just generate all windows
  // upfront to avoid the overhead of batched emission for smaller runs.
  static constexpr usize BATCH_THRESHOLD = 2 * BATCH_SIZE;
  const bool use_batching = num_total_windows > BATCH_THRESHOLD;

  // Track all emitted windows for variant flushing (indexed by genome index)
  std::vector<core::WindowPtr> windows;
  windows.reserve(num_total_windows);

  absl::FixedArray<bool> done_windows(num_total_windows);
  done_windows.fill(false);

  const auto send_qptr = std::make_shared<AsyncWorker::InputQueue>(num_total_windows);
  const auto recv_qptr = std::make_shared<AsyncWorker::OutputQueue>(num_total_windows);
  const moodycamel::ProducerToken producer_token(*send_qptr);

  usize batch_offset = 0;

  // Generates the next batch of windows from the builder, enqueues them for
  // workers, and appends them to the tracking vector. Callers are responsible
  // for checking whether more windows are needed before invoking.
  const auto feed_next_batch = [&window_builder, &batch_offset, &send_qptr, &producer_token, &windows]() {
    auto next_batch = window_builder.BuildWindowsBatch(batch_offset);
    if (!next_batch.empty()) {
      send_qptr->enqueue_bulk(producer_token, next_batch.begin(), next_batch.size());
      windows.insert(windows.end(), next_batch.begin(), next_batch.end());
    }
  };

  if (use_batching) {
    // Seed the queue with the first batch of windows
    feed_next_batch();
  } else {
    // Small run: generate all windows upfront, no batching overhead
    windows = window_builder.BuildWindows();
    batch_offset = num_total_windows;  // Mark all as emitted
    send_qptr->enqueue_bulk(producer_token, windows.begin(), windows.size());
  }

  // Launch worker threads
  std::vector<std::jthread> worker_threads;
  worker_threads.reserve(mParamsPtr->mNumWorkerThreads);
  const auto varstore = std::make_shared<core::VariantStore>();
  const auto vb_params = std::make_shared<const core::VariantBuilder::Params>(mParamsPtr->mVariantBuilder);
  const auto window_length = mParamsPtr->mWindowBuilder.mWindowLength;
  for (usize idx = 0; idx < mParamsPtr->mNumWorkerThreads; ++idx) {
    worker_threads.emplace_back(PipelineWorker, &producer_token, send_qptr, recv_qptr, varstore, vb_params, window_length);
  }

  static const auto all_windows_upto_idx_done = [](const absl::FixedArray<bool> &dw, const usize window_idx) -> bool {
    const auto *last_itr = window_idx >= dw.size() ? dw.cend() : dw.cbegin() + window_idx;
    return std::all_of(dw.cbegin(), last_itr, [](const bool is_done) { return is_done; });
  };

  const auto percent_done = [&num_total_windows](const usize ndone) -> f64 {
    return 100.0 * (static_cast<f64>(ndone) / static_cast<f64>(num_total_windows));
  };

  usize idx_to_flush = 0;
  usize num_completed = 0;
  core::AsyncWorker::Result async_worker_result;
  moodycamel::ConsumerToken result_consumer_token(*recv_qptr);

  auto stats = InitWindowStats();
  constexpr usize nbuffer_windows = 100;
  EtaTimer eta_timer(num_total_windows);

  // ---------------------------------------------------------------------------
  // Main pipeline loop: process results and dynamically feed new window batches
  // ---------------------------------------------------------------------------
  while (num_completed != num_total_windows) {
    if (!recv_qptr->try_dequeue(result_consumer_token, async_worker_result)) {
      using namespace std::chrono_literals;
      std::this_thread::sleep_for(1s);
      if (use_batching && batch_offset < num_total_windows && send_qptr->size_approx() < BATCH_SIZE) {
        feed_next_batch();
      }
      continue;
    }

    num_completed++;
    stats.at(async_worker_result.mStatus) += 1;
    done_windows[async_worker_result.mGenomeIdx] = true;
    const core::WindowPtr &curr_win = windows[async_worker_result.mGenomeIdx];
    const auto win_name = curr_win->ToSamtoolsRegion();
    const auto win_status = core::ToString(async_worker_result.mStatus);

    eta_timer.Increment();
    const auto elapsed_rt = absl::FormatDuration(absl::Trunc(timer.Runtime(), absl::Seconds(1)));
    const auto rem_rt = absl::FormatDuration(absl::Trunc(eta_timer.EstimatedEta(), absl::Seconds(1)));
    const auto win_rt = absl::FormatDuration(absl::Trunc(async_worker_result.mRuntime, absl::Microseconds(100)));

    LOG_INFO("Progress: {:>8.4f}% | Elapsed: {} | ETA: {} @ {:.2f}/s | {} done with {} in {}",
             percent_done(num_completed), elapsed_rt, rem_rt, eta_timer.RatePerSecond(), win_name, win_status, win_rt)

    if (all_windows_upto_idx_done(done_windows, idx_to_flush + nbuffer_windows)) {
      varstore->FlushVariantsBeforeWindow(*windows[idx_to_flush], output_vcf);
      idx_to_flush++;
    }

    if (use_batching && batch_offset < num_total_windows && send_qptr->size_approx() < BATCH_SIZE) {
      feed_next_batch();
    }
  }

  std::ranges::for_each(worker_threads, std::mem_fn(&std::jthread::request_stop));
  std::ranges::for_each(worker_threads, std::mem_fn(&std::jthread::join));

#ifdef LANCET_PROFILE_MODE
  ProfilerStop();
  ProfilerFlush();
#endif

  varstore->FlushAllVariantsInStore(output_vcf);
  output_vcf.Close();

  LogWindowStats(stats);
  const auto total_runtime = absl::FormatDuration(absl::Trunc(timer.Runtime(), absl::Milliseconds(1)));
  LOG_INFO("Successfully completed processing {} windows | Runtime={}", num_total_windows, total_runtime)
  std::exit(EXIT_SUCCESS);
}

// ---------------------------------------------------------------------------
// BuildVcfHeader
// ---------------------------------------------------------------------------

auto PipelineRunner::BuildVcfHeader(const CliParams &params) -> std::string {
  using namespace std::string_view_literals;
  // clang-format off
  static constexpr auto fstr_hdr = R"raw(##fileformat=VCFv4.3
##fileDate={RUN_TIMESTAMP}
##source=Lancet_{FULL_VERSION_TAG}
##commandLine="{FULL_COMMAND_USED}"
##reference="{REFERENCE_PATH}"
{CONTIG_HDR_LINES}##INFO=<ID=SHARED,Number=0,Type=Flag,Description="Variant ALT seen in both tumor & normal sample(s)">
##INFO=<ID=NORMAL,Number=0,Type=Flag,Description="Variant ALT seen only in normal samples(s)">
##INFO=<ID=TUMOR,Number=0,Type=Flag,Description="Variant ALT seen only in tumor sample(s)">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type. Possible values are SNV, INS, DEL and MNP">
##INFO=<ID=LENGTH,Number=1,Type=Integer,Description="Variant length in base pairs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype called at the variant site">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth per allele (REF, ALT1, ALT2, ...)">
##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Forward strand read depth per allele">
##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Reverse strand read depth per allele">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the variant site">
##FORMAT=<ID=RMQ,Number=R,Type=Float,Description="RMS mapping quality per allele">
##FORMAT=<ID=PBQ,Number=R,Type=Float,Description="Posterior base quality per allele (Bayesian aggregation)">
##FORMAT=<ID=SB,Number=R,Type=Float,Description="Strand bias ratio per allele (fwd/total)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality (second-lowest PL, capped at 99)">
)raw"sv;
  // clang-format on

  std::string contig_hdr_lines;
  static constexpr usize CONTIGS_BUFFER_SIZE = 524288;
  contig_hdr_lines.reserve(CONTIGS_BUFFER_SIZE);
  const hts::Reference ref(params.mVariantBuilder.mRdCollParams.mRefPath);
  for (const auto &chrom : ref.ListChroms()) {
    absl::StrAppend(&contig_hdr_lines, fmt::format("##contig=<ID={},length={}>\n", chrom.Name(), chrom.Length()));
  }

  auto full_hdr = fmt::format(
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      fstr_hdr, fmt::arg("RUN_TIMESTAMP", absl::FormatTime(absl::RFC3339_sec, absl::Now(), absl::LocalTimeZone())),
      fmt::arg("FULL_VERSION_TAG", LancetFullVersion()), fmt::arg("FULL_COMMAND_USED", params.mFullCmdLine),
      fmt::arg("REFERENCE_PATH", params.mVariantBuilder.mRdCollParams.mRefPath.string()),
      fmt::arg("CONTIG_HDR_LINES", contig_hdr_lines));

  const auto rc_sample_list = core::ReadCollector::BuildSampleNameList(params.mVariantBuilder.mRdCollParams);
  absl::StrAppend(&full_hdr, fmt::format("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
                                         absl::StrJoin(rc_sample_list, "\t")));

  return full_hdr;
}

// ---------------------------------------------------------------------------
// ValidateAndPopulateParams
// ---------------------------------------------------------------------------

void PipelineRunner::ValidateAndPopulateParams() {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mParamsPtr->mVariantBuilder.mSkipActiveRegion) return;

  using lancet::hts::Alignment;
  static constexpr usize NUM_READS_TO_PEEK = 1000;
  static const std::vector<std::string> tags{"MD"};
  const lancet::hts::Reference ref(mParamsPtr->mVariantBuilder.mRdCollParams.mRefPath);

  const auto is_md_missing = [&ref](const std::filesystem::path &aln_path) -> bool {
    usize peeked_read_count = 0;
    lancet::hts::Extractor extractor(aln_path, ref, Alignment::Fields::AUX_RGAUX, tags, true);
    for (const auto &aln : extractor) {
      // NOLINTBEGIN(readability-braces-around-statements)
      if (peeked_read_count > NUM_READS_TO_PEEK) break;
      if (aln.HasTag("MD")) return false;
      // NOLINTEND(readability-braces-around-statements)
      peeked_read_count++;
    }
    return true;
  };

  if (std::ranges::any_of(mParamsPtr->mVariantBuilder.mRdCollParams.mNormalPaths, is_md_missing)) {
    LOG_WARN("MD tag missing in normal BAM/CRAM. Turning off active region detection")
    mParamsPtr->mVariantBuilder.mSkipActiveRegion = true;
    return;
  }

  if (std::ranges::any_of(mParamsPtr->mVariantBuilder.mRdCollParams.mTumorPaths, is_md_missing)) {
    LOG_WARN("MD tag missing in tumor BAM/CRAM. Turning off active region detection")
    mParamsPtr->mVariantBuilder.mSkipActiveRegion = true;
    return;
  }
}

}  // namespace lancet::cli
