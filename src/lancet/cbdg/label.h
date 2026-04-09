#ifndef SRC_LANCET_CBDG_LABEL_H_
#define SRC_LANCET_CBDG_LABEL_H_

#include "lancet/base/types.h"

namespace lancet::cbdg {

class Label {
 public:
  // NOLINTNEXTLINE(cppcoreguidelines-use-enum-class)
  enum Tag : u8 { REFERENCE = 1, NORMAL = 2, TUMOR = 4 };

  Label() = default;
  constexpr explicit Label(Tag val) noexcept : mData(static_cast<u8>(val)) {}
  constexpr explicit Label(u8 data) noexcept : mData(data) {}

  void Merge(Tag const tag) { mData |= static_cast<u8>(tag); }
  void Merge(Label const& other) { mData |= other.mData; }

  [[nodiscard]] constexpr auto HasTag(Tag const tag) const noexcept -> bool {
    return (mData & static_cast<u8>(tag)) != 0;
  }

  [[nodiscard]] constexpr auto GetData() const noexcept -> u8 { return mData; }

  auto operator==(Label const& rhs) const noexcept -> bool { return mData == rhs.mData; }
  auto operator!=(Label const& rhs) const noexcept -> bool { return !(*this == rhs); }

  auto operator==(Tag const tag) const noexcept -> bool { return mData == static_cast<u8>(tag); }
  auto operator!=(Tag const tag) const noexcept -> bool { return !(*this == tag); }

 private:
  u8 mData = 0;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_LABEL_H_
