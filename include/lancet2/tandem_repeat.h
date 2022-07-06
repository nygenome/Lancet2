#pragma once

#include <string>
#include <string_view>

#include "lancet2/sized_ints.h"

namespace lancet2 {
struct TandemRepeatParams {
  u32 maxSTRUnitLength = 4;
  u32 minSTRUnitsToReport = 3;
  u32 minSTRLengthToReport = 7;
  u32 distFromSTR = 1;
};

struct TandemRepeatResult {
  bool foundSTR = false;
  usize strLength = 0;
  std::string strMotif = "";
};

auto FindTandemRepeat(std::string_view seq, usize pos, const TandemRepeatParams& params) -> TandemRepeatResult;
}  // namespace lancet2
