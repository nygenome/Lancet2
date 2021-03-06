#include "lancet2/run_pipeline.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <future>
#include <utility>
#include <vector>

#include "absl/container/fixed_array.h"
#include "lancet2/assert_macro.h"
#include "lancet2/fasta_reader.h"
#include "lancet2/hts_reader.h"
#include "lancet2/log_macros.h"
#include "lancet2/micro_assembler.h"
#include "lancet2/timer.h"
#include "lancet2/variant_store.h"
#include "lancet2/window_builder.h"
#include "spdlog/spdlog.h"

namespace lancet2 {
static inline auto GetSampleNames(const CliParams& p) -> std::vector<std::string> {
  HtsReader rdrN(p.normalPath, p.referencePath);
  HtsReader rdrT(p.tumorPath, p.referencePath);

  return {rdrN.GetSampleName(), rdrT.GetSampleName()};
}

static inline auto GetContigIDs(const CliParams& p) -> absl::flat_hash_map<std::string, i64> {
  return FastaReader(p.referencePath).GetContigIndexMap();
}

static inline auto RequiredBufferWindows(const CliParams& p) -> usize {
  // Number of windows ahead of current window to be done in order to flush current window
  const auto maxFlankLen = static_cast<double>(std::max(p.maxIndelLength, p.windowLength));
  const auto windowStep = static_cast<double>(WindowBuilder::StepSize(p.pctOverlap, p.windowLength));
  return static_cast<usize>(4.0 * std::ceil(maxFlankLen / windowStep));
}

void RunPipeline(std::shared_ptr<CliParams> params) {  // NOLINT
  Timer T;
  LOG_INFO("Starting main thread for processing lancet2 pipeline");
  if (!params->ValidateParams()) std::exit(EXIT_FAILURE);
  LOG_INFO("Successfully validated input command line parameters");

  if (!params->outGraphsDir.empty()) {
    std::filesystem::remove_all(params->outGraphsDir);
    std::filesystem::create_directory(params->outGraphsDir);
  }

  std::ofstream outVcf(params->outVcfPath, std::ios_base::out | std::ios_base::trunc);
  const auto vcfHdr = VariantStore::GetHeader(GetSampleNames(*params), *params);
  outVcf.write(vcfHdr.c_str(), vcfHdr.length());

  const auto contigIDs = GetContigIDs(*params);
  const auto allwindows = BuildWindows(contigIDs, *params);
  const auto numThreads = static_cast<usize>(params->numWorkerThreads);
  const auto paramsPtr = std::make_shared<const CliParams>(*params);
  const auto numBufWindows = RequiredBufferWindows(*paramsPtr);
  const auto vDBPtr = std::make_shared<VariantStore>(paramsPtr);

  LOG_INFO("Processing {} windows in {} microassembler thread(s)", allwindows.size(), params->numWorkerThreads);
  std::vector<std::future<void>> assemblers;
  assemblers.reserve(numThreads);

  const auto resultQueuePtr = std::make_shared<OutResultQueue>(allwindows.size());
  const auto windowQueuePtr = std::make_shared<InWindowQueue>(allwindows.size());
  moodycamel::ProducerToken windowProducerToken(*windowQueuePtr);
  windowQueuePtr->enqueue_bulk(windowProducerToken, allwindows.begin(), allwindows.size());

  static absl::FixedArray<bool> doneWindows(allwindows.size());
  doneWindows.fill(false);
  std::set_terminate([]() -> void {
    LOG_CRITICAL("Caught unexpected program termination call! Exiting abnormally...");
    for (bool& doneWindow : doneWindows) doneWindow = true;
    std::abort();
  });

  for (usize idx = 0; idx < numThreads; ++idx) {
    assemblers.emplace_back(std::async(
        std::launch::async, [&vDBPtr](std::unique_ptr<MicroAssembler> m) -> void { m->Process(vDBPtr); },
        std::make_unique<MicroAssembler>(windowQueuePtr, resultQueuePtr, paramsPtr)));
  }

  const auto anyWindowsRunning = []() -> bool {
    return std::any_of(doneWindows.cbegin(), doneWindows.cend(), [](const bool& wdone) { return !wdone; });
  };

  const auto allWindowsUptoDone = [](const usize win_idx) -> bool {
    const auto* lastItr = win_idx >= doneWindows.size() ? doneWindows.cend() : doneWindows.cbegin() + win_idx;
    return std::all_of(doneWindows.cbegin(), lastItr, [](const bool& wdone) { return wdone; });
  };

  usize idxToFlush = 0;
  usize numDone = 0;
  const auto numTotal = allwindows.size();
  const auto pctDone = [&numTotal](const usize done) -> double {
    return 100.0 * (static_cast<double>(done) / static_cast<double>(numTotal));
  };

  WindowResult result;
  moodycamel::ConsumerToken resultConsumerToken(*resultQueuePtr);
  while (anyWindowsRunning()) {
    resultQueuePtr->wait_dequeue(resultConsumerToken, result);

    numDone++;
    doneWindows[result.windowIdx] = true;
    const auto windowID = allwindows[result.windowIdx]->ToRegionString();
    LOG_INFO("Progress: {:>7.3f}% | {} processed in {}", pctDone(numDone), windowID, Humanized(result.runtime));

    if (allWindowsUptoDone(idxToFlush + numBufWindows)) {
      const auto flushed = vDBPtr->FlushWindow(*allwindows[idxToFlush], outVcf, contigIDs);
      if (flushed) {
        LOG_DEBUG("Flushed variants from {} to output vcf", allwindows[idxToFlush]->ToRegionString());
        outVcf.flush();
      }

      idxToFlush++;
    }
  }

  vDBPtr->FlushAll(outVcf, contigIDs);
  outVcf.close();

  // just to make sure futures get collected and threads released
  std::for_each(assemblers.begin(), assemblers.end(), [](std::future<void>& fut) { return fut.get(); });
  LOG_INFO("Successfully completed lancet2 pipeline | Runtime={}", T.HumanRuntime());
  std::exit(EXIT_SUCCESS);
}
}  // namespace lancet2
