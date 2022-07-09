#include "lancet2/barcode_set.h"

#include <algorithm>

namespace lancet2 {
void BarcodeSet::Merge(const BarcodeSet &other) {
  for (usize idx = 0; idx < 4; idx++) {
    const auto &otherItem = other.data.at(idx);
    data.at(idx).insert(otherItem.cbegin(), otherItem.cend());
  }
}

auto BarcodeSet::IsEmpty() const -> bool {
  return std::all_of(data.cbegin(), data.cend(), [](const auto &c) { return c.empty(); });
}

auto BarcodeSet::GetSize() const -> usize {
  usize result = 0;
  std::for_each(data.cbegin(), data.cend(), [&result](const auto &c) { result += c.size(); });
  return result;
}

auto BarcodeSet::AddBX(SampleLabel label, Strand s, std::string_view bx) -> bool {
  return data.at(ToIdx(label, s)).insert(std::string(bx)).second;
}

auto BarcodeSet::IsBXMissing(SampleLabel label, Strand s, std::string_view bx) const -> bool {
  const auto &c = data.at(ToIdx(label, s));
  return c.find(bx) == c.end();
}

auto BarcodeSet::IsBXMissing(SampleLabel label, std::string_view bx) const -> bool {
  return IsBXMissing(label, Strand::FWD, bx) || IsBXMissing(label, Strand::REV, bx);
}

auto BarcodeSet::BXCount(SampleLabel label, Strand s) const -> u16 {
  return static_cast<u16>(data.at(ToIdx(label, s)).size());
}

auto BarcodeSet::ToIdx(SampleLabel label, Strand s) -> usize {
  switch (label) {
    case SampleLabel::NORMAL:
      return s == Strand::FWD ? 0 : 1;
    case SampleLabel::TUMOR:
    default:
      return s == Strand::FWD ? 2 : 3;
  }
}

}  // namespace lancet2
