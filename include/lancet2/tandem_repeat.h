#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace lancet2 {
struct TandemRepeatParams {
  std::uint32_t maxSTRUnitLength = 4;
  std::uint32_t minSTRUnitsToReport = 3;
  std::uint32_t minSTRLengthToReport = 7;
  std::uint32_t distFromSTR = 1;
};

struct TandemRepeatResult {
  bool foundSTR = false;
  std::size_t strLength = 0;
  std::string strMotif = "";
};

auto FindTandemRepeat(std::string_view seq, std::size_t pos, const TandemRepeatParams& params) -> TandemRepeatResult;
}  // namespace lancet2
