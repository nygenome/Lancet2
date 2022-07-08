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

  // add unique bases first and then merge labels for k-1 overlapping bases
  usize buddyIdx = 0;
  if (dir == BuddyPosition::FRONT) {
    if (reverse_buddy) {
      bases.insert(bases.end(), buddy.bases.crbegin() + k - 1, buddy.bases.crend());

      // merge overlapping base labels
      buddyIdx = buddy.bases.size() - 1;
      for (usize idx = bases.size() - k; idx < bases.size(); ++idx) {
        bases.at(idx) = (bases.at(idx) | buddy.bases.at(buddyIdx));
        buddyIdx--;
      }
    } else {
      bases.insert(bases.end(), buddy.bases.cbegin() + k - 1, buddy.bases.cend());

      // merge overlapping base labels
      buddyIdx = 0;
      for (usize idx = bases.size() - k; idx < bases.size(); ++idx) {
        bases.at(idx) = (bases.at(idx) | buddy.bases.at(buddyIdx));
        buddyIdx++;
      }
    }
  } else {
    if (reverse_buddy) {
      bases.insert(bases.begin(), buddy.bases.crbegin(), buddy.bases.crend() - k + 1);

      // merge overlapping base labels
      buddyIdx = k - 1;
      for (usize idx = buddy.bases.size() - k; idx < buddy.bases.size() - 1; ++idx) {
        bases.at(idx) = (bases.at(idx) | buddy.bases.at(buddyIdx));
        buddyIdx--;
      }
    } else {
      bases.insert(bases.begin(), buddy.bases.cbegin(), buddy.bases.cend() - k + 1);

      // merge overlapping base labels
      buddyIdx = buddy.bases.size() - k;
      for (usize idx = buddy.bases.size() - k; idx < buddy.bases.size() - 1; ++idx) {
        bases.at(idx) = (bases.at(idx) | buddy.bases.at(buddyIdx));
        buddyIdx++;
      }
    }
  }
}

void NodeLabel::Push(KmerLabel label) {
  const auto setLabel = [&label](BaseLabel& base) { base.SetLabel(label, true); };
  std::for_each(bases.begin(), bases.end(), setLabel);
}

auto NodeLabel::UniqueLabelRatio(KmerLabel label) const -> double {
  const auto hasUniqLabel = [&label](const BaseLabel& base) { return base.IsLabelOnly(label); };
  const auto count = std::count_if(bases.cbegin(), bases.cend(), hasUniqLabel);
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

  if (hasRef && hasTmr && hasNml) return "lightblue";
  if (hasTmr && !hasNml) return "orangered";
  if (hasNml && !hasTmr) return hasRef ? "lightblue" : "royalblue";
  return "lightblue";
}
}  // namespace lancet2
