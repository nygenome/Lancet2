#include "lancet2/node_qual.h"

#include "lancet2/assert_macro.h"

namespace lancet2 {
NodeQual::NodeQual(usize count) { data.resize(count); }

void NodeQual::Push(absl::string_view sv) {
  LANCET_ASSERT(data.size() == sv.size());  // NOLINT

  for (usize idx = 0; idx < sv.size(); idx++) {
    data[idx].Add(static_cast<double>(sv[idx]));
  }
}

auto NodeQual::LowQualPositions(double max_bq) const -> std::vector<bool> {
  std::vector<bool> result;
  for (const auto& baseQual : data) result.emplace_back(baseQual.GetMean() < max_bq);
  return result;
}

auto NodeQual::HighQualPositions(double min_bq) const -> std::vector<bool> {
  std::vector<bool> result;
  for (const auto& baseQual : data) result.emplace_back(baseQual.GetMean() >= min_bq);
  return result;
}
}  // namespace lancet2
