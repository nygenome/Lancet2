#include "lancet/cli_params.h"

#include "lancet/contig_info.h"
#include "lancet/fasta_reader.h"
#include "lancet/hts_reader.h"
#include "spdlog/spdlog.h"

namespace lancet {
static const auto TagPresent = [](const CliParams& p, const char* tag) -> bool {
  return HasTag(p.tumorPath, p.referencePath, tag) || HasTag(p.normalPath, p.referencePath, tag);
};

static const auto RefContigsMatch = [](const CliParams& p) -> bool {
  FastaReader refFa(p.referencePath);
  HtsReader rdrT(p.tumorPath, p.referencePath);
  HtsReader rdrN(p.normalPath, p.referencePath);

  return CheckContigsMatch(rdrT.ContigsInfo(), refFa.ContigsInfo()) &&
         CheckContigsMatch(rdrN.ContigsInfo(), refFa.ContigsInfo());
};

auto CliParams::ValidateParams() -> bool {
  // ensure MD tag is present when active region is not turned off
  if (!activeRegionOff && !TagPresent(*this, "MD")) {
    SPDLOG_WARN("MD tag is missing from tumor and normal BAMs/CRAMs. Turning off active region detection.");
    activeRegionOff = true;
  }

  // ensure HP and BX tags are present when tenxMode is turned on
  if (tenxMode) {
    if (!TagPresent(*this, "BX")) {
      SPDLOG_WARN("BX tag is missing from tumor and normal BAMs/CRAMs. Turning off tenx-mode.");
      tenxMode = false;
    }

    if (!TagPresent(*this, "HP")) {
      SPDLOG_WARN("HP tag is missing from tumor and normal BAMs/CRAMs. Turning off tenx-mode.");
      tenxMode = false;
    }
  }

  if (!noCtgCheck && !RefContigsMatch(*this)) {
    SPDLOG_ERROR("Reference contigs in tumor/normal BAM/CRAMs do not match with reference FASTA.");
    return false;
  }

  return true;
}
}  // namespace lancet
