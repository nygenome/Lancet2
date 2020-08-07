#include "lancet/micro_assembler.h"

#include <algorithm>
#include <utility>

#include "absl/strings/str_format.h"
#include "lancet/fasta_reader.h"
#include "lancet/graph_builder.h"
#include "lancet/logger.h"
#include "lancet/read_extractor.h"
#include "lancet/utils.h"

namespace lancet {
void MicroAssembler::Process(const std::shared_ptr<VariantStore>& store) const {
  Timer T;
  FastaReader refRdr(params->referencePath);
  auto window = std::make_shared<RefWindow>();
  moodycamel::ConsumerToken ctok(*windowQPtr);

  while (windowQPtr->try_dequeue(ctok, window)) {
    T.Reset();
    const auto winIdx = window->WindowIndex();
    const auto regStr = window->ToRegionString();
    const auto regResult = refRdr.RegionSequence(window->ToGenomicRegion());

    if (!regResult.ok() && absl::IsFailedPrecondition(regResult.status()) && params->skipTruncSeq) {
      DebugLog("Skipping window %s with truncated reference sequence in fasta", regStr);
      resultNotifiers[winIdx]->SetResult(winIdx, T.Runtime());
      continue;
    }

    if (!regResult.ok()) {
      resultNotifiers[winIdx]->SetResult(winIdx, T.Runtime());
      WarnLog("Error processing %s: %s", regStr, regResult.status().message());
      continue;
    }

    window->SetSequence(regResult.ValueOrDie());
    const auto procStatus = ProcessWindow(std::const_pointer_cast<const RefWindow>(window), store);

    if (!procStatus.ok()) {
      const auto errMsg = absl::StrFormat("Error processing %s: %s", regStr, procStatus.message());
      WarnLog("%s", errMsg);
      resultNotifiers[winIdx]->SetResult(winIdx, T.Runtime());
      continue;
    }

    resultNotifiers[winIdx]->SetResult(winIdx, T.Runtime());
  }
}

auto MicroAssembler::ProcessWindow(const std::shared_ptr<const RefWindow>& w,
                                   const std::shared_ptr<VariantStore>& store) const -> absl::Status {
  const auto regionStr = w->ToRegionString();
  DebugLog("Starting to process %s in MicroAssembler", regionStr);
  if (ShouldSkipWindow(w)) return absl::OkStatus();

  ReadExtractor re(params, w->ToGenomicRegion());
  if (!params->activeRegionOff && !re.IsActiveRegion()) {
    DebugLog("Skipping %s since no evidence of mutation is found", regionStr);
    return absl::OkStatus();
  }

  const auto reads = re.Extract();
  GraphBuilder gb(w, absl::MakeConstSpan(reads), re.AverageCoverage(), params);
  auto graph = gb.BuildGraph(params->minKmerSize, params->maxKmerSize);
  graph->ProcessGraph({gb.RefData(SampleLabel::NORMAL), gb.RefData(SampleLabel::TUMOR)}, store);

  while (graph->ShouldIncrementK()) {
    if (gb.CurrentKmerSize() == params->maxKmerSize) {
      WarnLog("Skipping %s after trying max kmer=%d", regionStr, params->maxKmerSize);
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
    WarnLog("Skipping %s since it has only N bases in reference", regionStr);
    return true;
  }

  if (utils::HasRepeatKmer(refseq, params->maxKmerSize)) {
    WarnLog("Skipping %s since reference has repeat %d-mers", regionStr, params->maxKmerSize);
    return true;
  }

  return false;
}
}  // namespace lancet
