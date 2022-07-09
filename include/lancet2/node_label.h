#pragma once

#include <vector>

#include "lancet2/base_label.h"
#include "lancet2/core_enums.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class NodeLabel {
 public:
  explicit NodeLabel(usize node_len);
  NodeLabel() = delete;

  void MergeBuddy(const NodeLabel& buddy, BuddyPosition dir, bool reverse_buddy, usize k);

  void Push(KmerLabel label);

  [[nodiscard]] auto LabelRatio(KmerLabel label) const -> double;
  [[nodiscard]] auto HasLabel(KmerLabel label) const -> bool;
  [[nodiscard]] auto IsLabelOnly(KmerLabel label) const -> bool;

  [[nodiscard]] auto GetLength() const noexcept -> usize { return bases.size(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return bases.empty(); }

  [[nodiscard]] auto FillColor() const -> std::string;

  void Reserve(const usize count) { bases.reserve(count); }

 private:
  std::vector<BaseLabel> bases;
};
}  // namespace lancet2
