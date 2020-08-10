#include "lancet/node_qual.h"

#include "absl/types/span.h"
#include "lancet/assert_macro.h"
#include "lancet/merge_node_info.h"

namespace lancet {
NodeQual::NodeQual(std::size_t count) { data.resize(count); }

void NodeQual::MergeBuddy(const NodeQual& buddy, BuddyPosition dir, bool reverse_buddy, std::size_t k) {
  MergeNodeInfo(&data, absl::MakeConstSpan(buddy.data), dir, reverse_buddy, k);
}

void NodeQual::Push(absl::string_view sv) {
  LANCET_ASSERT(data.size() == sv.size());  // NOLINT

  for (std::size_t idx = 0; idx < sv.size(); idx++) {
    data[idx].Add(static_cast<double>(sv[idx]));
  }
}

auto NodeQual::LowQualPositions(double max_bq) const -> std::vector<bool> {
  std::vector<bool> result;
  for (const auto& baseQual : data) result.emplace_back(baseQual.Mean() < max_bq);
  return result;
}

auto NodeQual::HighQualPositions(double min_bq) const -> std::vector<bool> {
  std::vector<bool> result;
  for (const auto& baseQual : data) result.emplace_back(baseQual.Mean() >= min_bq);
  return result;
}
}  // namespace lancet
