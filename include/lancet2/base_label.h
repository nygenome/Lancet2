#pragma once

#include <bitset>
#include <initializer_list>
#include <string>

#include "lancet2/core_enums.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class BaseLabel {
 public:
  BaseLabel() = default;
  explicit BaseLabel(KmerLabel label);
  BaseLabel(std::initializer_list<KmerLabel> labels);

  [[nodiscard]] auto GetCount() const noexcept -> usize;

  void SetLabel(KmerLabel label, bool value = true);

  [[nodiscard]] auto HasLabel(KmerLabel label) const -> bool;
  [[nodiscard]] auto IsLabelOnly(KmerLabel label) const -> bool;

  auto operator&(const BaseLabel& other) const -> BaseLabel;
  auto operator&=(const BaseLabel& other) -> BaseLabel&;

  auto operator|(const BaseLabel& other) const -> BaseLabel;
  auto operator|=(const BaseLabel& other) -> BaseLabel&;

  auto operator^(const BaseLabel& other) const -> BaseLabel;
  auto operator^=(const BaseLabel& other) -> BaseLabel&;

  [[nodiscard]] auto ToString() const -> std::string;

 private:
  std::bitset<3> bits;

  [[nodiscard]] static auto ToBitPosition(KmerLabel label) -> usize;
};
}  // namespace lancet2
