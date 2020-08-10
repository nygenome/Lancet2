#pragma once

#include <cstddef>
#include <vector>

#include "absl/strings/string_view.h"
#include "lancet/core_enums.h"
#include "lancet/online_stats.h"

namespace lancet {
class NodeQual {
 public:
  explicit NodeQual(std::size_t count);
  NodeQual() = delete;

  void MergeBuddy(const NodeQual &buddy, BuddyPosition dir, bool reverse_buddy, std::size_t k);

  void Push(absl::string_view sv);

  [[nodiscard]] auto LowQualPositions(double max_bq) const -> std::vector<bool>;
  [[nodiscard]] auto HighQualPositions(double min_bq) const -> std::vector<bool>;

  [[nodiscard]] auto Length() const noexcept -> std::size_t { return data.size(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return data.empty(); }

  void Reserve(const std::size_t count) { data.reserve(count); }

 private:
  std::vector<OnlineStats> data;
};
}  // namespace lancet
