#ifndef SRC_LANCET_CBDG_LABEL_H_
#define SRC_LANCET_CBDG_LABEL_H_

#include "lancet/base/types.h"

namespace lancet::cbdg {

class Label {
 public:
  enum Tag : u8 { REFERENCE = 1, NORMAL = 2, TUMOR = 4 };

  Label() = default;
  constexpr Label(Tag val) noexcept : mData(static_cast<u8>(val)) {}
  constexpr Label(u8 data) noexcept : mData(data) {}

  void Merge(const Tag tag) { mData |= static_cast<u8>(tag); }
  void Merge(const Label& other) { mData |= other.mData; }

  [[nodiscard]] constexpr auto HasTag(const Tag tag) const noexcept -> bool {
    return (mData & static_cast<u8>(tag)) != 0;
  }

  [[nodiscard]] constexpr auto GetData() const noexcept -> u8 { return mData; }

  auto operator==(const Label& rhs) const noexcept -> bool { return mData == rhs.mData; }
  auto operator!=(const Label& rhs) const noexcept -> bool { return !(*this == rhs); }

  auto operator==(const Tag tag) const noexcept -> bool { return mData == static_cast<u8>(tag); }
  auto operator!=(const Tag tag) const noexcept -> bool { return !(*this == tag); }

 private:
  u8 mData = 0;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_LABEL_H_
