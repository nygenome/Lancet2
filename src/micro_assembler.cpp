#include "lancet/micro_assembler.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "absl/strings/str_format.h"
#include "lancet/graph_builder.h"
#include "lancet/logger.h"
#include "lancet/read_extractor.h"
#include "lancet/utils.h"

namespace lancet {
auto MicroAssembler::Process(const std::shared_ptr<VariantStore>& store) const -> absl::Status {
  DebugLog("Starting MicroAssembler to process %d windows", windows.size());

  for (std::size_t idx = 0; idx < windows.size(); ++idx) {
    const auto winIdx = windows[idx]->WindowIndex();
    const auto winStatus = ProcessWindow(windows[idx], store);

    if (!winStatus.ok()) {
      const auto regionStr = windows[idx]->ToRegionString();
      const auto errMsg = absl::StrFormat("Error processing %s in MicroAssembler", regionStr);
      FatalLog("%s: %s", errMsg, winStatus.message());
      throw std::runtime_error(winStatus.ToString());
    }

    notifiers[idx]->SetResult(winIdx);
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
