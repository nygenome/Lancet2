#include "lancet2/cli_params.h"

#include "lancet2/contig_info.h"
#include "lancet2/fasta_reader.h"
#include "lancet2/hts_reader.h"
#include "lancet2/log_macros.h"
#include "spdlog/spdlog.h"

namespace lancet2 {
static const auto TagPresent = [](const CliParams& p, const char* tag) -> bool {
  return TagPeekCheck(p.tumorPath, p.referencePath, tag) || TagPeekCheck(p.normalPath, p.referencePath, tag);
};

static const auto RefContigsMatch = [](const CliParams& p) -> bool {
  FastaReader refFa(p.referencePath);
  HtsReader rdrT(p.tumorPath, p.referencePath);
  HtsReader rdrN(p.normalPath, p.referencePath);

  return CheckContigsMatch(rdrT.GetContigs(), refFa.GetContigs()) &&
         CheckContigsMatch(rdrN.GetContigs(), refFa.GetContigs());
};

auto CliParams::ValidateParams() -> bool {
  // ensure MD tag is present when active region is not turned off
  if (!activeRegionOff && !TagPresent(*this, "MD")) {
    LOG_WARN("MD tag is missing from tumor and normal BAMs/CRAMs. Turning off active region detection.");
    activeRegionOff = true;
  }

  // ensure HP and BX tags are present when tenxMode is turned on
  if (tenxMode) {
    if (!TagPresent(*this, "BX")) {
      LOG_WARN("BX tag is missing from tumor and normal BAMs/CRAMs. Turning off tenx-mode.");
      tenxMode = false;
    }

    if (!TagPresent(*this, "HP")) {
      LOG_WARN("HP tag is missing from tumor and normal BAMs/CRAMs. Turning off tenx-mode.");
      tenxMode = false;
    }
  }

  if (!noCtgCheck && !RefContigsMatch(*this)) {
    LOG_ERROR("Reference contigs in tumor/normal BAM/CRAMs do not match with reference FASTA.");
    return false;
  }

  return true;
}
}  // namespace lancet2
