#include "lancet2/node_label.h"

#include <algorithm>

#include "absl/types/span.h"
#include "lancet2/merge_node_info.h"

namespace lancet2 {
NodeLabel::NodeLabel(usize node_len) { bases.resize(node_len); }

void NodeLabel::MergeBuddy(const NodeLabel& buddy, BuddyPosition dir, bool reverse_buddy, usize k) {
  if (bases.empty()) {
    reverse_buddy ? bases.insert(bases.end(), buddy.bases.crbegin(), buddy.bases.crend())
                  : bases.insert(bases.end(), buddy.bases.cbegin(), buddy.bases.cend());
    return;
  }

  MergeNodeInfo(&bases, absl::MakeConstSpan(buddy.bases), dir, reverse_buddy, k);
}

void NodeLabel::Push(KmerLabel label) {
  const auto setLabel = [&label](BaseLabel& base) { base.SetLabel(label, true); };
  std::for_each(bases.begin(), bases.end(), setLabel);
}

auto NodeLabel::LabelRatio(KmerLabel label) const -> double {
  const auto hasLabel = [&label](const BaseLabel& base) { return base.HasLabel(label); };
  const auto count = std::count_if(bases.cbegin(), bases.cend(), hasLabel);
  return static_cast<double>(count) / static_cast<double>(bases.size());
}

auto NodeLabel::HasLabel(KmerLabel label) const -> bool {
  const auto hasLabel = [&label](const BaseLabel& base) { return base.HasLabel(label); };
  return std::any_of(bases.cbegin(), bases.cend(), hasLabel);
}

auto NodeLabel::IsLabelOnly(KmerLabel label) const -> bool {
  const auto hasLabel = [&label](const BaseLabel& base) { return base.HasLabel(label); };
  return std::all_of(bases.cbegin(), bases.cend(), hasLabel);
}

auto NodeLabel::FillColor() const -> std::string {
  const auto hasRef = HasLabel(KmerLabel::REFERENCE);
  const auto hasTmr = HasLabel(KmerLabel::TUMOR);
  const auto hasNml = HasLabel(KmerLabel::NORMAL);

  if (hasRef || (hasTmr && hasNml)) return "royalblue"; // SHARED (or) REF
  if (hasTmr && !hasNml && !hasRef) return "orangered"; // TUMOR only
  if (hasNml && !hasTmr && !hasRef) return "limegreen"; // NORMAL only
  return "lightblue"; // Should not get here. But catch call if we get here.
}
}  // namespace lancet2
