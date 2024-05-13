#include "lancet/cli/pipeline_runner.h"

#include <algorithm>
#include <chrono>  // NOLINT(misc-include-cleaner)
#include <cmath>
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

// #ifndef LANCET_DEVELOP_MODE
// #include "gperftools/profiler.h"
// #endif

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
                    AsyncWorker::VariantStorePtr vstore, AsyncWorker::BuilderParamsPtr params) {
  // NOLINTEND(bugprone-easily-swappable-parameters,performance-unnecessary-value-param)
  // #ifndef LANCET_DEVELOP_MODE
  // NOLINTNEXTLINE(readability-braces-around-statements)
  // if (ProfilingIsEnabledForAllThreads() != 0) ProfilerRegisterThread();
  // #endif
  auto worker =
      std::make_unique<AsyncWorker>(std::move(in_queue), std::move(out_queue), std::move(vstore), std::move(params));
  worker->Process(std::move(stop_token), *in_token);
}

namespace lancet::cli {

PipelineRunner::PipelineRunner(std::shared_ptr<CliParams> params) : mParamsPtr(std::move(params)) {
  // #ifndef LANCET_DEVELOP_MODE
  //   setenv("CPUPROFILE_PER_THREAD_TIMERS", "1", 1);
  //   setenv("CPUPROFILE_FREQUENCY", "10000", 1);
  //   const auto timestamp = absl::FormatTime("%Y%m%d%ET%H%M%S", absl::Now(), absl::LocalTimeZone());
  //   const auto fname = fmt::format("Lancet.cpu_profile.{}.bin", timestamp);
  //   ProfilerStart(fname.c_str());
  // #endif
}

void PipelineRunner::Run() {
  Timer timer;
  static thread_local const auto tid = std::this_thread::get_id();
  LOG_INFO("Using main thread {:#x} to synchronize variant calling pipeline", absl::Hash<std::thread::id>()(tid))

  ValidateAndPopulateParams();
  if (!mParamsPtr->mVariantBuilder.mOutGraphsDir.empty()) {
    // Set out graphs directory parameter as well and create new out graphs root diretory
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
  const auto windows = BuildWindows(*mParamsPtr);
  LOG_INFO("Processing {} window(s) with {} VariantBuilder thread(s)", windows.size(), mParamsPtr->mNumWorkerThreads)

  const auto num_total_windows = windows.size();
  static absl::FixedArray<bool> done_windows(num_total_windows);
  done_windows.fill(false);

  const auto send_qptr = std::make_shared<AsyncWorker::InputQueue>(windows.size());
  const auto recv_qptr = std::make_shared<AsyncWorker::OutputQueue>(windows.size());
  const moodycamel::ProducerToken producer_token(*send_qptr);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overflow"
  send_qptr->enqueue_bulk(producer_token, windows.begin(), windows.size());
#pragma GCC diagnostic pop

  std::vector<std::jthread> worker_threads;
  worker_threads.reserve(mParamsPtr->mNumWorkerThreads);
  const auto varstore = std::make_shared<core::VariantStore>();
  const auto vb_params = std::make_shared<const core::VariantBuilder::Params>(mParamsPtr->mVariantBuilder);
  for (usize idx = 0; idx < mParamsPtr->mNumWorkerThreads; ++idx) {
    worker_threads.emplace_back(PipelineWorker, &producer_token, send_qptr, recv_qptr, varstore, vb_params);
  }

  static const auto all_windows_upto_idx_done = [](const usize window_idx) -> bool {
    const auto *last_itr = window_idx >= done_windows.size() ? done_windows.cend() : done_windows.cbegin() + window_idx;
    return std::all_of(done_windows.cbegin(), last_itr, [](const bool is_window_done) { return is_window_done; });
  };

  static const auto percent_done = [&num_total_windows](const usize ndone) -> f64 {
    return 100.0 * (static_cast<f64>(ndone) / static_cast<f64>(num_total_windows));
  };

  usize idx_to_flush = 0;
  usize num_completed = 0;
  core::AsyncWorker::Result async_worker_result;
  moodycamel::ConsumerToken result_consumer_token(*recv_qptr);

  auto stats = InitWindowStats();
  constexpr usize nbuffer_windows = 100;
  EtaTimer eta_timer(num_total_windows);

  while (num_completed != num_total_windows) {
    if (!recv_qptr->try_dequeue(result_consumer_token, async_worker_result)) {
      using namespace std::chrono_literals;
      std::this_thread::sleep_for(1s);
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

    if (all_windows_upto_idx_done(idx_to_flush + nbuffer_windows)) {
      varstore->FlushVariantsBeforeWindow(*windows[idx_to_flush], output_vcf);
      idx_to_flush++;
    }
  }

  std::ranges::for_each(worker_threads, std::mem_fn(&std::jthread::request_stop));
  std::ranges::for_each(worker_threads, std::mem_fn(&std::jthread::join));

  // #ifndef LANCET_DEVELOP_MODE
  //   ProfilerStop();
  //   ProfilerFlush();
  // #endif

  varstore->FlushAllVariantsInStore(output_vcf);
  output_vcf.Close();

  LogWindowStats(stats);
  const auto total_runtime = absl::FormatDuration(absl::Trunc(timer.Runtime(), absl::Milliseconds(1)));
  LOG_INFO("Successfully completed processing {} windows | Runtime={}", num_total_windows, total_runtime)
  std::exit(EXIT_SUCCESS);
}

auto PipelineRunner::BuildWindows(const CliParams &params) -> std::vector<core::WindowPtr> {
  core::WindowBuilder window_builder(params.mVariantBuilder.mRdCollParams.mRefPath, params.mWindowBuilder);
  window_builder.AddBatchRegions(absl::MakeConstSpan(params.mInRegions));
  window_builder.AddBatchRegions(params.mBedFile);

  if (window_builder.IsEmpty()) {
    LOG_WARN("No input regions provided to build windows. Using contigs in reference as input regions")
    window_builder.AddAllReferenceRegions();
  }

  return window_builder.BuildWindows();
}

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
##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant ALT seen near an identified STR site">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type. Possible values are SNV, INS, DEL and MNP">
##INFO=<ID=LENGTH,Number=1,Type=Integer,Description="Variant length in base pairs">
##INFO=<ID=KMERLEN,Number=1,Type=Integer,Description="K-mer length used to assemble the locus">
##INFO=<ID=STR_LEN,Number=1,Type=Integer,Description="If variant ALT is near STR, lists length of the STR unit">
##INFO=<ID=STR_MOTIF,Number=1,Type=String,Description="If variant ALT is near STR, lists motif of the STR unit">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype called at the variant site">
##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Number of reads supporting REF and ALT alleles">
##FORMAT=<ID=ADF,Number=2,Type=Integer,Description="Number of reads supporting REF and ALT alleles on forward strand">
##FORMAT=<ID=ADR,Number=2,Type=Integer,Description="Number of reads supporting REF and ALT alleles on reverse strand">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Read depth in the sample at the variant site">
##FORMAT=<ID=WDC,Number=1,Type=Float,Description="Window read depth after downsampling and read filters">
##FORMAT=<ID=WTC,Number=1,Type=Float,Description="Window read depth before downsampling and read filters">
##FORMAT=<ID=PRF,Number=1,Type=Float,Description="Fraction of reads in the window that pass read quality filters">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="ALT allele frequency in the sample at the variant site">
##FORMAT=<ID=RAQS,Number=4,Type=Integer,Description="REF allele quality stats - Min, Median, Max, MAD">
##FORMAT=<ID=AAQS,Number=4,Type=Integer,Description="ALT allele quality stats - Min, Median, Max, MAD">
##FORMAT=<ID=RMQS,Number=4,Type=Integer,Description="REF mapping quality stats - Min, Median, Max, MAD">
##FORMAT=<ID=AMQS,Number=4,Type=Integer,Description="ALT mapping quality stats - Min, Median, Max, MAD">
##FORMAT=<ID=RAPDS,Number=4,Type=Integer,Description="REF aln scores pct difference stats - Min, Median, Max, MAD">
##FORMAT=<ID=AAPDS,Number=4,Type=Integer,Description="ALT aln scores pct difference stats - Min, Median, Max, MAD">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred-scaled genotype quality for the sample">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized phred-scaled likelihoods for all genotypes">
)raw"sv;
  // clang-format on

  static const auto should_exclude_chrom = [](std::string_view chrom) -> bool {
    return chrom == "MT" || chrom == "chrM" || absl::StartsWith(chrom, "GL") || absl::StartsWith(chrom, "chrUn") ||
           absl::StartsWith(chrom, "chrEBV") || absl::StartsWith(chrom, "HLA-") || absl::EndsWith(chrom, "_random") ||
           absl::EndsWith(chrom, "_alt") || absl::EndsWith(chrom, "_decoy");
  };

  std::string contig_hdr_lines;
  static constexpr usize CONTIGS_BUFFER_SIZE = 524288;
  contig_hdr_lines.reserve(CONTIGS_BUFFER_SIZE);
  const hts::Reference ref(params.mVariantBuilder.mRdCollParams.mRefPath);
  for (const auto &chrom : ref.ListChroms()) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (should_exclude_chrom(chrom.Name())) continue;
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
