#include "lancet2/base_label.h"

namespace lancet2 {
BaseLabel::BaseLabel(KmerLabel label) { SetLabel(label); }

BaseLabel::BaseLabel(std::initializer_list<KmerLabel> labels) {
  for (const auto& label : labels) SetLabel(label);
}

auto BaseLabel::ToBitPosition(KmerLabel label) -> usize { return static_cast<usize>(static_cast<u8>(label)); }

auto BaseLabel::GetCount() const noexcept -> usize { return bits.count(); }

void BaseLabel::SetLabel(KmerLabel label, bool value) { bits.set(ToBitPosition(label), value); }

auto BaseLabel::HasLabel(KmerLabel label) const -> bool { return bits.test(ToBitPosition(label)); }
auto BaseLabel::IsLabelOnly(KmerLabel label) const -> bool { return GetCount() == 1 && HasLabel(label); }

auto BaseLabel::operator&(const BaseLabel& other) const -> BaseLabel {
  BaseLabel result;
  result.bits = this->bits & other.bits;
  return result;
}

auto BaseLabel::operator&=(const BaseLabel& other) -> BaseLabel& {
  bits &= other.bits;
  return *this;
}

auto BaseLabel::operator|(const BaseLabel& other) const -> BaseLabel {
  BaseLabel result;
  result.bits = this->bits | other.bits;
  return result;
}

auto BaseLabel::operator|=(const BaseLabel& other) -> BaseLabel& {
  bits |= other.bits;
  return *this;
}

auto BaseLabel::operator^(const BaseLabel& other) const -> BaseLabel {
  BaseLabel result;
  result.bits = this->bits ^ other.bits;
  return result;
}

auto BaseLabel::operator^=(const BaseLabel& other) -> BaseLabel& {
  bits ^= other.bits;
  return *this;
}

auto BaseLabel::ToString() const -> std::string {
  std::string result;
  if (HasLabel(KmerLabel::REFERENCE)) result += "R";
  if (HasLabel((KmerLabel::NORMAL))) result += "N";
  if (HasLabel(KmerLabel::TUMOR)) result += "T";
  return result;
}
}  // namespace lancet2
