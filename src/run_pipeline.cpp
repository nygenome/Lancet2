#include "lancet/run_pipeline.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <future>
#include <vector>

#include "absl/container/fixed_array.h"
#include "absl/status/status.h"
#include "lancet/fasta_reader.h"
#include "lancet/hts_reader.h"
#include "lancet/logger.h"
#include "lancet/micro_assembler.h"
#include "lancet/timer.h"
#include "lancet/utils.h"
#include "lancet/variant_store.h"
#include "lancet/vcf_writer.h"
#include "lancet/window_builder.h"

namespace lancet {
static inline auto GetSampleNames(const CliParams& p) -> std::vector<std::string> {
  HtsReader rdrN(p.normalPath, p.referencePath);
  const auto resultN = rdrN.SampleNames();
  assert(resultN.size() == 1);  // NOLINT

  HtsReader rdrT(p.tumorPath, p.referencePath);
  const auto resultT = rdrT.SampleNames();
  assert(resultT.size() == 1);  // NOLINT

  if (resultN.size() != 1 || resultT.size() != 1) {
    throw std::runtime_error("expected both tumor and normal BAM/CRAMs to have one sample name");
  }

  return {resultN[0], resultT[0]};
}

static inline auto GetContigIDs(const CliParams& p) -> absl::flat_hash_map<std::string, std::int64_t> {
  return FastaReader(p.referencePath).ContigIDs();
}

void RunPipeline(std::shared_ptr<CliParams> params) {  // NOLINT
  Timer T;
  InfoLog("Starting main thread for processing lancet pipeline");
  if (!params->ValidateParams()) std::exit(EXIT_FAILURE);
  InfoLog("Successfully validated input command line parameters");

  if (!params->outGraphsDir.empty()) {
    const auto result = utils::MakeDir(params->outGraphsDir);
    if (!result.ok()) {
      FatalLog("Could not create output graphs dir: %s; %s", params->outGraphsDir, result.message());
      std::exit(EXIT_FAILURE);
    }
  }

  VcfWriter outVcf(params->outVcfPath);
  if (!outVcf.Write(VariantStore::BuildVcfHeader(GetSampleNames(*params), *params)).ok()) {
    FatalLog("Could not write header to output vcf: %s", params->outVcfPath);
    std::exit(EXIT_FAILURE);
  }

  const auto contigIDs = GetContigIDs(*params);
  const auto allwindows = BuildWindows(contigIDs, *params);
  InfoLog("Processing %d windows in %d worker thread(s)", allwindows.size(), params->numWorkerThreads);

  std::vector<std::shared_ptr<Notification<std::size_t>>> notifiers(allwindows.size());
  for (std::size_t idx = 0; idx < allwindows.size(); ++idx) {
    notifiers[idx] = std::make_shared<Notification<std::size_t>>();
  }

  auto perTWindows = utils::RoundRobinSplit(absl::MakeConstSpan(allwindows), params->numWorkerThreads);
  auto perTNotifiers = utils::RoundRobinSplit(absl::MakeConstSpan(notifiers), params->numWorkerThreads);

  const auto paramsPtr = std::make_shared<const CliParams>(*params);
  const auto storePtr = std::make_shared<VariantStore>(allwindows.size(), paramsPtr);

  std::vector<std::future<absl::Status>> assemblers;
  assemblers.reserve(perTWindows.size());

  for (std::size_t idx = 0; idx < perTWindows.size(); ++idx) {
    auto threadWindows = absl::MakeSpan(perTWindows[idx]);
    auto threadNotifiers = absl::MakeSpan(perTNotifiers[idx]);

    assemblers.emplace_back(std::async(
        std::launch::async,
        [&storePtr](std::unique_ptr<MicroAssembler> m) -> absl::Status { return m->Process(storePtr); },
        std::make_unique<MicroAssembler>(threadWindows, threadNotifiers, paramsPtr)));
  }

  const auto anyWindowRunning = [&notifiers]() -> bool {
    return std::any_of(notifiers.cbegin(), notifiers.cend(),
                       [](const MicroAssembler::NotificationPtr& n) { return n->IsNotDone(); });
  };

  const auto allWindowsUptoDone = [&notifiers](const std::size_t win_idx) -> bool {
    auto lastItr = win_idx >= notifiers.size() ? notifiers.cend() : notifiers.cbegin() + win_idx;
    return std::all_of(notifiers.cbegin(), lastItr,
                       [](const MicroAssembler::NotificationPtr& n) { return n->IsDone(); });
  };

  const auto maxFlankLen = static_cast<double>(std::max(params->maxIndelLength, params->windowLength));
  const auto windowStep = static_cast<double>(WindowBuilder::StepSize(params->pctOverlap, params->windowLength));
  const auto numBufferWindows = static_cast<std::size_t>(3.0 * std::ceil(maxFlankLen / windowStep));

  absl::FixedArray<bool> alreadyLogged(allwindows.size());
  alreadyLogged.fill(false);

  std::size_t idxToFlush = 0;
  std::size_t numDone = 0;
  const auto numTotal = allwindows.size();
  const auto pctDone = [&numTotal](const std::size_t done) -> double {
    return 100.0 * (static_cast<double>(done) / static_cast<double>(numTotal));
  };

  while (anyWindowRunning()) {
    if (idxToFlush == allwindows.size()) break;

    for (const auto& notifier : notifiers) {
      if (notifier->IsNotDone() || alreadyLogged[notifier->WaitForResult()]) continue;

      numDone++;
      const auto winIdx = notifier->WaitForResult();
      alreadyLogged[winIdx] = true;

      const auto windowID = allwindows[winIdx]->ToRegionString();
      InfoLog("Progress: %0.3f%% | Done processing %s | %d of %d completed", pctDone(numDone), windowID, numDone,
              numTotal);
    }

    if (allWindowsUptoDone(idxToFlush + numBufferWindows) && storePtr->FlushWindow(idxToFlush, &outVcf, contigIDs)) {
      DebugLog("Flushed variants from %s to output vcf", allwindows[idxToFlush]->ToRegionString());
      idxToFlush++;
    }
  }

  storePtr->FlushAll(&outVcf, contigIDs);
  outVcf.Close();

  // just to make sure futures get collected and threads released
  for (auto& assembler : assemblers) {
    const auto status = assembler.get();
    if (!status.ok()) {
      FatalLog("Exited MicroAssembler thread with error: %s", status.message());
      std::exit(EXIT_FAILURE);
    }
  }

  InfoLog("Successfully completed lancet pipeline | Runtime=%s", T.HumanRuntime());
  std::exit(EXIT_SUCCESS);
}
}  // namespace lancet
