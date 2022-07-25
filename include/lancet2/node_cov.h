#pragma once

#include "lancet2/base_cov.h"
#include "lancet2/core_enums.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class NodeCov {
 public:
  NodeCov() = default;

  void MergeBuddy(const NodeCov& buddy, usize nodeLen, usize buddyLen, usize k);

  [[nodiscard]] auto StrandCov(SampleLabel label, Strand s) const -> u16;
  [[nodiscard]] auto TotalCov(SampleLabel label) const -> u16;

  void Update(SampleLabel label, Strand s);
  void Update(u16 val, SampleLabel label, Strand s);

  void Reset();

 private:
  u16 cntTumorFwd = 0;
  u16 cntTumorRev = 0;
  u16 cntNormalFwd = 0;
  u16 cntNormalRev = 0;
};
}  // namespace lancet2
