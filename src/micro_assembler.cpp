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
auto MicroAssembler::Process(const std::shared_ptr<VariantStore>& store) const -> absl::Status {
  DebugLog("Starting MicroAssembler to process %d windows", windows.size());

  FastaReader refRdr(params->referencePath);
  for (std::size_t idx = 0; idx < windows.size(); ++idx) {
    const auto regStr = windows[idx]->ToRegionString();
    const auto regResult = refRdr.RegionSequence(windows[idx]->ToGenomicRegion());

    if (!regResult.ok() && absl::IsFailedPrecondition(regResult.status()) && params->skipTruncSeq) {
      DebugLog("Skipping window %s with truncated reference sequence in fasta", regStr);
      continue;
    }

    if (!regResult.ok()) {
      notifiers[idx]->SetErrorMsg(regResult.status().ToString());
      return regResult.status();
    }

    windows[idx]->SetSequence(regResult.ValueOrDie());
    const auto winStatus = ProcessWindow(std::const_pointer_cast<const RefWindow>(windows[idx]), store);
    if (!winStatus.ok()) {
      const auto errMsg = absl::StrFormat("Error processing %s: %s", regStr, winStatus.message());
      notifiers[idx]->SetErrorMsg(errMsg);
      return absl::InternalError(errMsg);
    }

    notifiers[idx]->SetResult(windows[idx]->WindowIndex());
  }

  return absl::OkStatus();
}

auto MicroAssembler::ProcessWindow(const std::shared_ptr<const RefWindow>& w,
                                   const std::shared_ptr<VariantStore>& store) const -> absl::Status {
  Timer T;
  const auto regionStr = w->ToRegionString();
  DebugLog("Starting to process %s in MicroAssembler", regionStr);
  if (ShouldSkipWindow(w, &T)) return absl::OkStatus();

  ReadExtractor re(params, w->ToGenomicRegion());
  if (!params->activeRegionOff && !re.IsActiveRegion()) {
    DebugLog("Skipping %s since no evidence of mutation is found", regionStr);
    DebugLog("Done processing %s in MicroAssembler | Runtime=%s", regionStr, T.HumanRuntime());
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

  DebugLog("Done processing %s in MicroAssembler | Runtime=%s", regionStr, T.HumanRuntime());
  return absl::OkStatus();
}

auto MicroAssembler::ShouldSkipWindow(const std::shared_ptr<const RefWindow>& w, Timer* T) const -> bool {
  const auto refseq = w->SeqView();
  const auto regionStr = w->ToRegionString();

  if (static_cast<std::size_t>(std::count(refseq.begin(), refseq.end(), 'N')) == refseq.length()) {
    WarnLog("Skipping %s since it has only N bases in reference | Runtime=%s", regionStr, T->HumanRuntime());
    return true;
  }

  if (utils::HasRepeatKmer(refseq, params->maxKmerSize)) {
    WarnLog("Skipping %s since reference has repeat %d-mers | Runtime=%s", regionStr, params->maxKmerSize,
            T->HumanRuntime());
    return true;
  }

  return false;
}
}  // namespace lancet
