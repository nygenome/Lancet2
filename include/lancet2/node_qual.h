#pragma once

#include <vector>

#include "absl/strings/string_view.h"
#include "lancet2/core_enums.h"
#include "lancet2/online_stats.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class NodeQual {
 public:
  explicit NodeQual(usize count);
  NodeQual() = delete;

  void MergeBuddy(const NodeQual &buddy, BuddyPosition dir, bool reverse_buddy, usize k);

  void Push(absl::string_view sv);

  [[nodiscard]] auto LowQualPositions(double max_bq) const -> std::vector<bool>;
  [[nodiscard]] auto HighQualPositions(double min_bq) const -> std::vector<bool>;

  [[nodiscard]] auto GetLength() const noexcept -> usize { return data.size(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return data.empty(); }

  void Reserve(const usize count) { data.reserve(count); }

 private:
  std::vector<OnlineStats> data;
};
}  // namespace lancet2
