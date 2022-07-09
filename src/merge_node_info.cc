#include "lancet2/merge_node_info.h"

#include "lancet2/assert_macro.h"
#include "lancet2/utils.h"

namespace lancet2 {
auto CanMergeSeqs(std::string_view source, std::string_view buddy, BuddyPosition merge_dir, bool reverse_buddy, usize k)
    -> bool {
  LANCET_ASSERT(!source.empty());  // NOLINT
  const auto bSeq = reverse_buddy ? utils::RevComp(buddy) : std::string(buddy);
  std::string_view buddySeq = bSeq;
  if (merge_dir == BuddyPosition::BACK) {
    return source.substr(0, k - 1) == buddySeq.substr(buddySeq.length() - k + 1, k - 1);
  }

  return source.substr(source.length() - k + 1, k - 1) == buddySeq.substr(0, k - 1);
}

void MergeKmerSeqs(std::string* source, std::string_view buddy, BuddyPosition dir, bool reverse_buddy, usize k) {
  LANCET_ASSERT(source != nullptr && !source->empty());  // NOLINT
#ifndef NDEBUG
  LANCET_ASSERT(CanMergeSeqs(*source, buddy, dir, reverse_buddy, k));  // NOLINT
#endif

  switch (dir) {
    case BuddyPosition::FRONT:
      reverse_buddy ? utils::PushRevCompSeq(source, source->end(), buddy, k - 1, buddy.length())
                    : utils::PushSeq(source, source->end(), buddy, k - 1, buddy.length());
      break;
    case BuddyPosition::BACK:
    default:
      reverse_buddy ? utils::PushRevCompSeq(source, source->begin(), buddy, 0, buddy.length() - k + 1)
                    : utils::PushSeq(source, source->begin(), buddy, 0, buddy.length() - k + 1);
      break;
  }
}
}  // namespace lancet2
