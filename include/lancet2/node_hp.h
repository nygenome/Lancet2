#pragma once

#include <vector>

#include "lancet2/base_hp.h"
#include "lancet2/core_enums.h"
#include "lancet2/node_cov.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class NodeHP {
 public:
  explicit NodeHP(const NodeCov& node_cov);
  NodeHP() = default;

  void MergeBuddy(const NodeHP& buddy, BuddyPosition dir, bool reverse_buddy, usize k);

  void Update(usize hp, SampleLabel label, const std::vector<bool>& bq_pass);
  void Update(usize hp, SampleLabel label, usize base_position);

  [[nodiscard]] auto BaseHPs(SampleLabel label) const noexcept -> std::vector<BaseHP> {
    return label == SampleLabel::TUMOR ? tmrHPs : nmlHPs;
  }

  [[nodiscard]] auto At(SampleLabel label, usize pos) -> BaseHP& {
    return label == SampleLabel::TUMOR ? tmrHPs.at(pos) : nmlHPs.at(pos);
  }

  [[nodiscard]] auto At(SampleLabel label, usize pos) const -> const BaseHP& {
    return label == SampleLabel::TUMOR ? tmrHPs.at(pos) : nmlHPs.at(pos);
  }

  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return tmrHPs.empty(); }
  [[nodiscard]] auto GetSize() const noexcept -> usize { return tmrHPs.size(); }

  void Clear();
  void Reverse();

  void Reserve(const usize count) {
    tmrHPs.reserve(count);
    nmlHPs.reserve(count);
  }

 private:
  std::vector<BaseHP> tmrHPs;
  std::vector<BaseHP> nmlHPs;
};
}  // namespace lancet2
