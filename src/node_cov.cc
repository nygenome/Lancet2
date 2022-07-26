#include "lancet2/node_cov.h"

#include <cmath>

namespace lancet2 {
auto NodeCov::StrandCov(SampleLabel label, Strand s) const -> u16 {
  if (label == SampleLabel::TUMOR) return s == Strand::FWD ? cntTumorFwd : cntTumorRev;
  return s == Strand::FWD ? cntNormalFwd : cntNormalRev;
}

auto NodeCov::TotalCov(SampleLabel label) const -> u16 {
  return StrandCov(label, Strand::FWD) + StrandCov(label, Strand::REV);
}

void NodeCov::Update(SampleLabel label, Strand s) {
  if (label == SampleLabel::TUMOR) {
    s == Strand::FWD ? cntTumorFwd++ : cntTumorRev++;
  } else {
    s == Strand::FWD ? cntNormalFwd++ : cntNormalRev++;
  }
}

void NodeCov::Update(u16 val, SampleLabel label, Strand s) {
  if (label == SampleLabel::TUMOR) {
    s == Strand::FWD ? cntTumorFwd = val : cntTumorRev = val;
  } else {
    s == Strand::FWD ? cntNormalFwd = val : cntNormalRev = val;
  }
}

void NodeCov::Reset() {
  cntTumorFwd = 0;
  cntTumorRev = 0;
  cntNormalFwd = 0;
  cntNormalRev = 0;
}
}  // namespace lancet2
