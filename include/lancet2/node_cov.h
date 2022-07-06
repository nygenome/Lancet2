#pragma once

#include <vector>

#include "lancet2/base_cov.h"
#include "lancet2/core_enums.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class NodeCov {
 public:
  explicit NodeCov(usize count) : tmrBases(count), nmlBases(count) {}
  NodeCov() = default;

  void MergeBuddy(const NodeCov& buddy, BuddyPosition dir, bool reverse_buddy, usize k);

  [[nodiscard]] auto StrandCov(SampleLabel label, Strand s) const -> u16;
  [[nodiscard]] auto TotalCov(SampleLabel label) const -> u16;

  void Update(SampleLabel label, Strand s, usize pos);
  void Update(SampleLabel label, Strand s, const std::vector<bool>& bq_pass);
  void Update(u16 val, SampleLabel label, Strand s, const std::vector<bool>& bq_pass);

  [[nodiscard]] auto BaseCovs(SampleLabel label) const noexcept -> std::vector<BaseCov> {
    return label == SampleLabel::TUMOR ? tmrBases : nmlBases;
  }

  auto At(SampleLabel label, usize pos) -> BaseCov& {
    return label == SampleLabel::TUMOR ? tmrBases.at(pos) : nmlBases.at(pos);
  }

  [[nodiscard]] auto At(SampleLabel label, usize pos) const -> const BaseCov& {
    return label == SampleLabel::TUMOR ? tmrBases.at(pos) : nmlBases.at(pos);
  }

  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return tmrBases.empty(); }
  [[nodiscard]] auto Size() const noexcept -> usize { return tmrBases.size(); }

  void Clear();
  void Reverse();

  void Reserve(const usize count) {
    tmrBases.reserve(count);
    nmlBases.reserve(count);
  }

 private:
  u16 cntTumorFwd = 0;
  u16 cntTumorRev = 0;
  u16 cntNormalFwd = 0;
  u16 cntNormalRev = 0;

  std::vector<BaseCov> tmrBases;
  std::vector<BaseCov> nmlBases;
};
}  // namespace lancet2
