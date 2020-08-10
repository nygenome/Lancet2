#include "lancet/micro_assembler.h"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <thread>

#include "absl/hash/hash.h"
#include "absl/strings/str_format.h"
#include "lancet/fasta_reader.h"
#include "lancet/graph_builder.h"
#include "lancet/timer.h"
#include "lancet/utils.h"
#include "spdlog/spdlog.h"

namespace lancet {
void MicroAssembler::Process(const std::shared_ptr<VariantStore>& store) const {
  static thread_local const auto tid = std::this_thread::get_id();
  SPDLOG_INFO("Starting MicroAssembler thread {:#x}", absl::Hash<std::thread::id>()(tid));

  Timer T;
  FastaReader refRdr(params->referencePath);
  ReadExtractor readExtractor(params);
  auto window = std::make_shared<RefWindow>();
  moodycamel::ConsumerToken windowConsumerToken(*windowQPtr);
  moodycamel::ProducerToken resultProducerToken(*resultQPtr);
  std::size_t numProcessed = 0;

  while (windowQPtr->try_dequeue(windowConsumerToken, window)) {
    T.Reset();
    const auto winIdx = window->WindowIndex();
    const auto regStr = window->ToRegionString();
    const auto regResult = refRdr.RegionSequence(window->ToGenomicRegion());
    numProcessed++;

    if (!regResult.ok() && absl::IsFailedPrecondition(regResult.status())) {
      SPDLOG_DEBUG("Skipping window {} with truncated reference sequence in fasta", regStr);
      resultQPtr->enqueue(resultProducerToken, WindowResult{T.Runtime(), winIdx});
      continue;
    }

    if (!regResult.ok()) {
      SPDLOG_ERROR("Error processing window {}: {}", regStr, regResult.status().message());
      resultQPtr->enqueue(resultProducerToken, WindowResult{T.Runtime(), winIdx});
      continue;
    }

    try {
      window->SetSequence(regResult.ValueOrDie());
      const auto windowStatus = ProcessWindow(&readExtractor, std::const_pointer_cast<const RefWindow>(window), store);
      if (!windowStatus.ok()) SPDLOG_ERROR("Error processing window {}: {}", regStr, windowStatus.message());
    } catch (const std::exception& exception) {
      SPDLOG_ERROR("Error processing window {}: {}", regStr, exception.what());
    } catch (...) {
      SPDLOG_ERROR("Error processing window {}: unknown exception caught", regStr);
    }

    resultQPtr->enqueue(resultProducerToken, WindowResult{T.Runtime(), winIdx});
  }

  SPDLOG_INFO("Processed {} windows in MicroAssembler thread {:#x}", numProcessed, absl::Hash<std::thread::id>()(tid));
  SPDLOG_INFO("Exiting MicroAssembler thread {:#x}", absl::Hash<std::thread::id>()(tid));
}

auto MicroAssembler::ProcessWindow(ReadExtractor* re, const std::shared_ptr<const RefWindow>& w,
                                   const std::shared_ptr<VariantStore>& store) const -> absl::Status {
  const auto regionStr = w->ToRegionString();
  SPDLOG_DEBUG("Starting to process {} in MicroAssembler", regionStr);
  if (ShouldSkipWindow(w)) return absl::OkStatus();

  re->SetTargetRegion(w->ToGenomicRegion());
  if (!params->activeRegionOff && !re->IsActiveRegion()) {
    SPDLOG_DEBUG("Skipping {} since no evidence of mutation is found", regionStr);
    return absl::OkStatus();
  }

  const auto reads = re->Extract();
  GraphBuilder gb(w, absl::MakeConstSpan(reads), re->AverageCoverage(), params);
  auto graph = gb.BuildGraph(params->minKmerSize, params->maxKmerSize);
  graph->ProcessGraph({gb.RefData(SampleLabel::NORMAL), gb.RefData(SampleLabel::TUMOR)}, store);

  while (graph->ShouldIncrementK()) {
    if (gb.CurrentKmerSize() == params->maxKmerSize) {
      SPDLOG_DEBUG("Skipping {} after trying to build graph with max k={}", regionStr, params->maxKmerSize);
      break;
    }

    graph = gb.BuildGraph(gb.CurrentKmerSize() + 2, params->maxKmerSize);
    graph->ProcessGraph({gb.RefData(SampleLabel::NORMAL), gb.RefData(SampleLabel::TUMOR)}, store);
  }

  return absl::OkStatus();
}

auto MicroAssembler::ShouldSkipWindow(const std::shared_ptr<const RefWindow>& w) const -> bool {
  const auto refseq = w->SeqView();
  const auto regionStr = w->ToRegionString();

  if (static_cast<std::size_t>(std::count(refseq.begin(), refseq.end(), 'N')) == refseq.length()) {
    SPDLOG_DEBUG("Skipping {} since it has only N bases in reference", regionStr);
    return true;
  }

  if (utils::HasRepeatKmer(refseq, params->maxKmerSize)) {
    SPDLOG_DEBUG("Skipping {} since reference has repeat {}-mers", regionStr, params->maxKmerSize);
    return true;
  }

  return false;
}
}  // namespace lancet
