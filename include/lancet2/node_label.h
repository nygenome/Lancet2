#pragma once

#include <cstddef>
#include <vector>

#include "lancet2/base_label.h"
#include "lancet2/core_enums.h"

namespace lancet2 {
class NodeLabel {
 public:
  explicit NodeLabel(std::size_t node_len);
  NodeLabel() = delete;

  void MergeBuddy(const NodeLabel& buddy, BuddyPosition dir, bool reverse_buddy, std::size_t k);

  void Push(KmerLabel label);

  [[nodiscard]] auto LabelRatio(KmerLabel label) const -> double;
  [[nodiscard]] auto HasLabel(KmerLabel label) const -> bool;
  [[nodiscard]] auto IsLabelOnly(KmerLabel label) const -> bool;

  [[nodiscard]] auto Length() const noexcept -> std::size_t { return bases.size(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return bases.empty(); }

  [[nodiscard]] auto FillColor() const -> std::string;

  void Reserve(const std::size_t count) { bases.reserve(count); }

 private:
  std::vector<BaseLabel> bases;
};
}  // namespace lancet2
