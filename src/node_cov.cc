#include "lancet2/node_cov.h"

#include <cmath>

namespace lancet2 {
void NodeCov::MergeBuddy(const NodeCov &buddy, usize len, usize buddyLen, usize k) {
  const auto aUniqLen = static_cast<double>(len - k + 1);
  const auto bUniqLen = static_cast<double>(buddyLen - k + 1);
  const auto abCommonLen = static_cast<double>(k - 1);
  const auto totalLen = static_cast<double>(len + buddyLen - k + 1);

  // Weighted harmonic mean (per base) and round to u16
  cntTumorFwd =
      static_cast<u16>(std::round(totalLen / ((aUniqLen / static_cast<double>(cntTumorFwd)) +
                                              (abCommonLen / static_cast<double>(cntTumorFwd + buddy.cntTumorFwd)) +
                                              (bUniqLen / static_cast<double>(buddy.cntTumorFwd)))));
  cntTumorRev =
      static_cast<u16>(std::round(totalLen / ((aUniqLen / static_cast<double>(cntTumorRev)) +
                                              (abCommonLen / static_cast<double>(cntTumorRev + buddy.cntTumorRev)) +
                                              (bUniqLen / static_cast<double>(buddy.cntTumorRev)))));
  cntNormalFwd =
      static_cast<u16>(std::round(totalLen / ((aUniqLen / static_cast<double>(cntNormalFwd)) +
                                              (abCommonLen / static_cast<double>(cntNormalFwd + buddy.cntNormalFwd)) +
                                              (bUniqLen / static_cast<double>(buddy.cntNormalFwd)))));
  cntNormalRev =
      static_cast<u16>(std::round(totalLen / ((aUniqLen / static_cast<double>(cntNormalRev)) +
                                              (abCommonLen / static_cast<double>(cntNormalRev + buddy.cntNormalRev)) +
                                              (bUniqLen / static_cast<double>(buddy.cntNormalRev)))));
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
