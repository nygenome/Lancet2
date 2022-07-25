#include "lancet2/node_cov.h"

#include <algorithm>
#include <cmath>

#include "absl/types/span.h"
#include "lancet2/merge_node_info.h"

namespace lancet2 {
void NodeCov::MergeBuddy(const NodeCov& buddy, usize nodeLen, usize buddyLen, usize k) {
  const auto totalLen = nodeLen + buddyLen - k + 1;
  const auto c1Ratio = static_cast<double>(nodeLen) / static_cast<double>(totalLen);
  const auto c2Ratio = static_cast<double>(buddyLen) / static_cast<double>(totalLen);

  cntTumorFwd = std::ceil(((cntTumorFwd * c1Ratio) + (buddy.cntTumorFwd * c2Ratio)));     // NOLINT
  cntTumorRev = std::ceil(((cntTumorRev * c1Ratio) + (buddy.cntTumorRev * c2Ratio)));     // NOLINT
  cntNormalFwd = std::ceil(((cntNormalFwd * c1Ratio) + (buddy.cntNormalFwd * c2Ratio)));  // NOLINT
  cntNormalRev = std::ceil(((cntNormalRev * c1Ratio) + (buddy.cntNormalRev * c2Ratio)));  // NOLINT
}

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
    s == Strand::FWD ? cntTumorFwd++ : cntTumorRev++;
  } else {
    s == Strand::FWD ? cntNormalFwd++ : cntNormalRev++;
  }
}

void NodeCov::Reset() {
  cntTumorFwd = 0;
  cntTumorRev = 0;
  cntNormalFwd = 0;
  cntNormalRev = 0;
}
}  // namespace lancet2
