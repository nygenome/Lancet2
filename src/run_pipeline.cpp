#include "lancet/run_pipeline.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <future>
#include <vector>

#include "absl/container/fixed_array.h"
#include "lancet/assert_macro.h"
#include "lancet/fasta_reader.h"
#include "lancet/hts_reader.h"
#include "lancet/micro_assembler.h"
#include "lancet/timer.h"
#include "lancet/utils.h"
#include "lancet/variant_store.h"
#include "lancet/window_builder.h"
#include "spdlog/spdlog.h"

namespace lancet {
static inline auto GetSampleNames(const CliParams& p) -> std::vector<std::string> {
  HtsReader rdrN(p.normalPath, p.referencePath);
  const auto resultN = rdrN.SampleNames();
  LANCET_ASSERT(resultN.size() == 1);  // NOLINT

  HtsReader rdrT(p.tumorPath, p.referencePath);
  const auto resultT = rdrT.SampleNames();
  LANCET_ASSERT(resultT.size() == 1);  // NOLINT

  if (resultN.size() != 1 || resultT.size() != 1) {
    throw std::runtime_error("expected both tumor and normal BAM/CRAMs to have one sample name");
  }

  return {resultN[0], resultT[0]};
}

static inline auto GetContigIDs(const CliParams& p) -> absl::flat_hash_map<std::string, std::int64_t> {
  return FastaReader(p.referencePath).ContigIDs();
}

static inline auto RequiredBufferWindows(const CliParams& p) -> std::size_t {
  // Number of windows ahead of current window to be done in order to flush current window
  const auto maxFlankLen = static_cast<double>(std::max(p.maxIndelLength, p.windowLength));
  const auto windowStep = static_cast<double>(WindowBuilder::StepSize(p.pctOverlap, p.windowLength));
  return static_cast<std::size_t>(3.0 * std::ceil(maxFlankLen / windowStep));
}

void RunPipeline(std::shared_ptr<CliParams> params) {  // NOLINT
  Timer T;
  SPDLOG_INFO("Starting main thread for processing lancet pipeline");
  if (!params->ValidateParams()) std::exit(EXIT_FAILURE);
  SPDLOG_INFO("Successfully validated input command line parameters");

  if (!params->outGraphsDir.empty()) {
    const auto result = utils::MakeDir(params->outGraphsDir);
    if (!result.ok()) {
      SPDLOG_CRITICAL("Could not create output graphs dir: {}; {}", params->outGraphsDir, result.message());
      std::exit(EXIT_FAILURE);
    }
  }

  std::ofstream outVcf(params->outVcfPath, std::ios_base::out | std::ios_base::trunc);
  const auto outHdr = VariantStore::BuildVcfHeader(GetSampleNames(*params), *params);
  outVcf.write(outHdr.c_str(), outHdr.length());

  const auto contigIDs = GetContigIDs(*params);
  const auto allwindows = BuildWindows(contigIDs, *params);
  const auto numThreads = static_cast<std::size_t>(params->numWorkerThreads);
  const auto paramsPtr = std::make_shared<const CliParams>(*params);
  const auto numBufWindows = RequiredBufferWindows(*paramsPtr);
  const auto vDBPtr = std::make_shared<VariantStore>(allwindows.size(), paramsPtr);

  SPDLOG_INFO("Processing {} windows in {} microassembler thread(s)", allwindows.size(), params->numWorkerThreads);
  std::vector<std::future<void>> assemblers;
  assemblers.reserve(numThreads);

  const auto resultQueuePtr = std::make_shared<OutResultQueue>(allwindows.size());
  const auto windowQueuePtr = std::make_shared<InWindowQueue>(allwindows.size());
  moodycamel::ProducerToken producerToken(*windowQueuePtr);
  windowQueuePtr->enqueue_bulk(producerToken, allwindows.begin(), allwindows.size());

  for (std::size_t idx = 0; idx < numThreads; ++idx) {
    assemblers.emplace_back(std::async(
        std::launch::async, [&vDBPtr](std::unique_ptr<MicroAssembler> m) -> void { return m->Process(vDBPtr); },
        std::make_unique<MicroAssembler>(windowQueuePtr, resultQueuePtr, paramsPtr)));
  }

  absl::FixedArray<bool> doneWindows(allwindows.size());
  doneWindows.fill(false);

  const auto anyWindowsRunning = [&doneWindows]() -> bool {
    return std::any_of(doneWindows.cbegin(), doneWindows.cend(), [](const bool& wdone) { return !wdone; });
  };

  const auto allWindowsUptoDone = [&doneWindows](const std::size_t win_idx) -> bool {
    const auto* lastItr = win_idx >= doneWindows.size() ? doneWindows.cend() : doneWindows.cbegin() + win_idx;
    return std::all_of(doneWindows.cbegin(), lastItr, [](const bool& wdone) { return wdone; });
  };

  std::size_t idxToFlush = 0;
  std::size_t numDone = 0;
  const auto numTotal = allwindows.size();
  const auto pctDone = [&numTotal](const std::size_t done) -> double {
    return 100.0 * (static_cast<double>(done) / static_cast<double>(numTotal));
  };

  WindowResult result;
  moodycamel::ConsumerToken consumerToken(*resultQueuePtr);
  while (anyWindowsRunning()) {
    resultQueuePtr->wait_dequeue(consumerToken, result);

    numDone++;
    doneWindows[result.windowIdx] = true;
    const auto windowID = allwindows[result.windowIdx]->ToRegionString();
    SPDLOG_INFO("Progress: {:03.3f}% | {} of {} done | Window {} processed in {}", pctDone(numDone), numDone, numTotal,
                windowID, Humanized(result.runtime));

    if (allWindowsUptoDone(idxToFlush + numBufWindows)) {
      const auto flushed = vDBPtr->FlushWindow(idxToFlush, outVcf, contigIDs);
      if (flushed) {
        SPDLOG_DEBUG("Flushed variants from {} to output vcf", allwindows[idxToFlush]->ToRegionString());
        outVcf.flush();
      }

      idxToFlush++;
    }
  }

  vDBPtr->FlushAll(outVcf, contigIDs);
  outVcf.close();

  // just to make sure futures get collected and threads released
  std::for_each(assemblers.begin(), assemblers.end(), [](std::future<void>& fut) { return fut.get(); });
  SPDLOG_INFO("Successfully completed lancet pipeline | Runtime={}", T.HumanRuntime());
  std::exit(EXIT_SUCCESS);
}
}  // namespace lancet
