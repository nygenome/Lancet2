#include "lancet2/micro_assembler.h"

#include <algorithm>
#include <exception>
#include <thread>

#include "absl/hash/hash.h"
#include "absl/strings/str_format.h"
#include "lancet2/fasta_reader.h"
#include "lancet2/graph_builder.h"
#include "lancet2/log_macros.h"
#include "lancet2/timer.h"
#include "lancet2/utils.h"
#include "spdlog/spdlog.h"

namespace lancet2 {
void MicroAssembler::Process(const std::shared_ptr<VariantStore>& store) {
  static thread_local const auto tid = std::this_thread::get_id();
  LOG_INFO("Started MicroAssembler thread {:#x}", absl::Hash<std::thread::id>()(tid));

  variants.reserve(2048);
  results.reserve(1024);

  Timer T;
  FastaReader refRdr(params->referencePath);
  ReadExtractor readExtractor(params);
  auto window = std::make_shared<RefWindow>();
  moodycamel::ConsumerToken windowConsumerToken(*windowQPtr);
  moodycamel::ProducerToken resultProducerToken(*resultQPtr);
  usize numProcessed = 0;

  while (windowQPtr->try_dequeue(windowConsumerToken, window)) {
    TryFlush(store, resultProducerToken);
    T.Reset();

    const auto winIdx = window->GetWindowIndex();
    const auto regStr = window->ToRegionString();
    const auto regResult = refRdr.GetRegionSeq(window->ToSamtoolsRegion());
    numProcessed++;

    if (!regResult.ok() && absl::IsFailedPrecondition(regResult.status())) {
      LOG_DEBUG("Skipping window {} with truncated reference sequence in fasta", regStr);
      results.emplace_back(WindowResult{T.GetRuntime(), winIdx});
      continue;
    }

    if (!regResult.ok()) {
      LOG_ERROR("Error processing window {}: {}", regStr, regResult.status().message());
      results.emplace_back(WindowResult{T.GetRuntime(), winIdx});
      continue;
    }

    try {
      window->SetSequence(regResult.value());
      const auto windowStatus = ProcessWindow(&readExtractor, std::const_pointer_cast<const RefWindow>(window));
      if (!windowStatus.ok()) LOG_ERROR("Error processing window {}: {}", regStr, windowStatus.message());
    } catch (const std::exception& exception) {
      LOG_ERROR("Error processing window {}: {}", regStr, exception.what());
    } catch (...) {
      LOG_ERROR("Error processing window {}: unknown exception caught", regStr);
    }

    results.emplace_back(WindowResult{T.GetRuntime(), winIdx});
  }

  ForceFlush(store, resultProducerToken);
  LOG_INFO("Done processing {} windows in MicroAssembler thread {:#x}", numProcessed,
           absl::Hash<std::thread::id>()(tid));
}

auto MicroAssembler::ProcessWindow(ReadExtractor* re, const std::shared_ptr<const RefWindow>& w) -> absl::Status {
  const auto regionStr = w->ToRegionString();
  LOG_DEBUG("Starting to process {} in MicroAssembler", regionStr);
  if (ShouldSkipWindow(w)) return absl::OkStatus();

  const auto scanResult = re->ScanRegion(w->ToSamtoolsRegion());
  if (!params->activeRegionOff && !scanResult.HasMutationEvidence) {
    LOG_DEBUG("Skipping {} since no evidence of mutation is found", regionStr);
    return absl::OkStatus();
  }

  const auto needToDownsample = scanResult.AverageCoverage <= params->maxWindowCov;
  const auto fraction = needToDownsample ? 1.0 : params->maxWindowCov / scanResult.AverageCoverage;
  const auto extractedCov = needToDownsample ? params->maxWindowCov : scanResult.AverageCoverage;
  const auto reads = re->ExtractReads(w->ToSamtoolsRegion(), fraction);

  GraphBuilder gb(w, absl::MakeConstSpan(reads), extractedCov, params);
  auto graph = gb.BuildGraph(params->minKmerSize, params->maxKmerSize);
  graph->ProcessGraph({gb.RefData(SampleLabel::NORMAL), gb.RefData(SampleLabel::TUMOR)}, &variants);

  while (graph->ShouldIncrementK()) {
    if (gb.CurrentKmerSize() == params->maxKmerSize) {
      LOG_DEBUG("Skipping {} after trying to build graph with max k={}", regionStr, params->maxKmerSize);
      break;
    }

    graph = gb.BuildGraph(gb.CurrentKmerSize() + 2, params->maxKmerSize);
    graph->ProcessGraph({gb.RefData(SampleLabel::NORMAL), gb.RefData(SampleLabel::TUMOR)}, &variants);
  }

  return absl::OkStatus();
}

auto MicroAssembler::ShouldSkipWindow(const std::shared_ptr<const RefWindow>& w) const -> bool {
  const auto refseq = w->GetSeqView();
  const auto regionStr = w->ToRegionString();

  if (static_cast<usize>(std::count(refseq.begin(), refseq.end(), 'N')) == refseq.length()) {
    LOG_DEBUG("Skipping {} since it has only N bases in reference", regionStr);
    return true;
  }

  if (utils::HasRepeatKmer(refseq, params->maxKmerSize)) {
    LOG_DEBUG("Skipping {} since reference has repeat {}-mers", regionStr, params->maxKmerSize);
    return true;
  }

  return false;
}

void MicroAssembler::TryFlush(const std::shared_ptr<VariantStore>& store, const moodycamel::ProducerToken& token) {
  const auto addedVariants = store->TryAddVariants(variants);
  if (addedVariants) {
    resultQPtr->enqueue_bulk(token, results.cbegin(), results.size());
    variants.clear();
    results.clear();
  }
}

void MicroAssembler::ForceFlush(const std::shared_ptr<VariantStore>& store, const moodycamel::ProducerToken& token) {
  store->ForceAddVariants(variants);
  resultQPtr->enqueue_bulk(token, results.cbegin(), results.size());
  variants.clear();
  results.clear();
}
}  // namespace lancet2
