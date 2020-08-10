#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "lancet/base_cov.h"
#include "lancet/core_enums.h"

namespace lancet {
class NodeCov {
 public:
  explicit NodeCov(std::size_t count) : tmrBases(count), nmlBases(count) {}
  NodeCov() = default;

  void MergeBuddy(const NodeCov& buddy, BuddyPosition dir, bool reverse_buddy, std::size_t k);

  [[nodiscard]] auto StrandCov(SampleLabel label, Strand s) const -> std::uint16_t;
  [[nodiscard]] auto TotalCov(SampleLabel label) const -> std::uint16_t;

  void Update(SampleLabel label, Strand s, std::size_t pos);
  void Update(SampleLabel label, Strand s, const std::vector<bool>& bq_pass);
  void Update(std::uint16_t val, SampleLabel label, Strand s, const std::vector<bool>& bq_pass);

  [[nodiscard]] auto BaseCovs(SampleLabel label) const noexcept -> std::vector<BaseCov> {
    return label == SampleLabel::TUMOR ? tmrBases : nmlBases;
  }

  auto At(SampleLabel label, std::size_t pos) -> BaseCov& {
    return label == SampleLabel::TUMOR ? tmrBases.at(pos) : nmlBases.at(pos);
  }

  [[nodiscard]] auto At(SampleLabel label, std::size_t pos) const -> const BaseCov& {
    return label == SampleLabel::TUMOR ? tmrBases.at(pos) : nmlBases.at(pos);
  }

  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return tmrBases.empty(); }
  [[nodiscard]] auto Size() const noexcept -> std::size_t { return tmrBases.size(); }

  void Clear();
  void Reverse();

  void Reserve(const std::size_t count) {
    tmrBases.reserve(count);
    nmlBases.reserve(count);
  }

 private:
  std::uint16_t cntTumorFwd = 0;
  std::uint16_t cntTumorRev = 0;
  std::uint16_t cntNormalFwd = 0;
  std::uint16_t cntNormalRev = 0;

  std::vector<BaseCov> tmrBases;
  std::vector<BaseCov> nmlBases;
};
}  // namespace lancet
