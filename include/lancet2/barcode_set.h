#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <string_view>

#include "absl/container/flat_hash_set.h"
#include "lancet2/core_enums.h"

namespace lancet2 {
class BarcodeSet {
 public:
  BarcodeSet() = default;

  void Merge(const BarcodeSet& other);

  [[nodiscard]] auto IsEmpty() const -> bool;
  [[nodiscard]] auto Size() const -> std::size_t;

  auto AddBX(SampleLabel label, Strand s, std::string_view bx) -> bool;  // NOLINT
  [[nodiscard]] auto IsBXMissing(SampleLabel label, Strand s, std::string_view bx) const -> bool;
  [[nodiscard]] auto IsBXMissing(SampleLabel label, std::string_view bx) const -> bool;
  [[nodiscard]] auto BXCount(SampleLabel label, Strand s) const -> std::uint16_t;

 private:
  using Container = absl::flat_hash_set<std::string>;
  std::array<Container, 4> data;

  [[nodiscard]] static auto ToIdx(SampleLabel label, Strand s) -> std::size_t;
};
}  // namespace lancet2
