#include "lancet2/bgzf_ostream.h"

namespace lancet2 {
extern "C" {
#include "htslib/tbx.h"
}

auto BgzfOstream::open(const std::string &path, BgzfFormat ofmt) -> bool {
  outFmt = ofmt;
  auto result = buf.open(path, "w");
  rdbuf(&buf);
  return result;
}

void BgzfOstream::close() {
  buf.close();
  if (outFmt != BgzfFormat::UNSPECIFIED) build_index();
}

void BgzfOstream::build_index() {
  switch (outFmt) {
    case BgzfFormat::VCF:
      tbx_index_build(buf.FileName.c_str(), 0, &tbx_conf_vcf);
      break;
    case BgzfFormat::GFF:
      tbx_index_build(buf.FileName.c_str(), 0, &tbx_conf_gff);
      break;
    case BgzfFormat::BED:
      tbx_index_build(buf.FileName.c_str(), 0, &tbx_conf_bed);
      break;
    default:
      return;
  }
}
}  // namespace lancet2
