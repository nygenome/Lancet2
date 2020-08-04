#pragma once

#include <cstddef>
#include <vector>

#include "lancet/base_hp.h"
#include "lancet/core_enums.h"
#include "lancet/node_cov.h"

namespace lancet {
class NodeHP {
 public:
  explicit NodeHP(const NodeCov& node_cov);
  NodeHP() = default;

  void MergeBuddy(const NodeHP& buddy, BuddyPosition dir, bool reverse_buddy, std::size_t k);

  void Update(std::size_t hp, SampleLabel label, const std::vector<std::size_t>& bq_pass);
  void Update(std::size_t hp, SampleLabel label, std::size_t base_position);

  [[nodiscard]] auto BaseHPs(SampleLabel label) const noexcept -> std::vector<BaseHP> {
    return label == SampleLabel::TUMOR ? tmrHPs : nmlHPs;
  }

  [[nodiscard]] auto At(SampleLabel label, std::size_t pos) -> BaseHP& {
    return label == SampleLabel::TUMOR ? tmrHPs.at(pos) : nmlHPs.at(pos);
  }

  [[nodiscard]] auto At(SampleLabel label, std::size_t pos) const -> const BaseHP& {
    return label == SampleLabel::TUMOR ? tmrHPs.at(pos) : nmlHPs.at(pos);
  }

  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return tmrHPs.empty(); }
  [[nodiscard]] auto Size() const noexcept -> std::size_t { return tmrHPs.size(); }

  void Clear();
  void Reverse();

  void Reserve(const std::size_t count) {
    tmrHPs.reserve(count);
    nmlHPs.reserve(count);
  }

 private:
  std::vector<BaseHP> tmrHPs;
  std::vector<BaseHP> nmlHPs;
};
}  // namespace lancet
