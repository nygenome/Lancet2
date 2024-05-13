#include "lancet/hts/iterator.h"

extern "C" {
#include "htslib/sam.h"
}

#include <stdexcept>

namespace lancet::hts {

auto Iterator::operator==(const Iterator& rhs) const -> bool { return mParsedAln == rhs.mParsedAln; }

auto Iterator::operator!=(const Iterator& rhs) const -> bool { return !(rhs == *this); }

auto Iterator::operator++() -> Iterator& {
  if (mRawFilePtr != nullptr && mRawHdrPtr != nullptr && mRawItrPtr != nullptr && mRawAlnPtr != nullptr) {
    FetchNextAlignment();
  }

  return *this;
}

auto Iterator::operator++(int) -> Iterator {
  auto old_val = *this;
  operator++();
  return old_val;
}

// NOLINTNEXTLINE(misc-no-recursion)
void Iterator::FetchNextAlignment() {
  const auto next_result = sam_itr_next(mRawFilePtr, mRawItrPtr, mRawAlnPtr);
  if (next_result < -1) {
    throw std::runtime_error("Could not fetch next alignment from BAM/CRAM iterator");
  }

  if (next_result == -1) {
    // No more data in hts_itr. We clear all fields in alignment so
    // that we match with end() Iterator during comparison
    mParsedAln.ClearAllFields();
    return;
  }

  if (PassesFilter()) {
    // Read passed all filters provided in the filter expression
    mParsedAln.PopulateRequestedFields(mRawAlnPtr, mFieldsNeeded, mTagsNeeded);
    return;
  }

  // Alignment did not pass filters, go fetch next alignment
  FetchNextAlignment();
}

auto Iterator::PassesFilter() const -> bool {
  if (mRawFiltrPtr == nullptr) return true;  // NOLINT(readability-braces-around-statements)

  // Filter expression present, so we need to apply filters to check
  // if the read passes filters before populating Alignment
  const auto check_result = sam_passes_filter(mRawHdrPtr, mRawAlnPtr, mRawFiltrPtr);
  if (check_result < 0) {
    // Error checking if alignment passed filters
    throw std::runtime_error("Could not apply filter expression to alignment");
  }

  return check_result == 1;
}

}  // namespace lancet::hts
