#ifndef SRC_LANCET_BASE_FIND_STR_H_
#define SRC_LANCET_BASE_FIND_STR_H_

#include <string>
#include <string_view>

#include "lancet/base/types.h"

struct StrParams {
  static constexpr usize DEFAULT_MAX_STR_UNIT_LENGTH = 4;
  static constexpr usize DEFAULT_MIN_STR_UNITS = 3;
  static constexpr usize DEFAULT_MIN_STR_LENGTH = 7;
  static constexpr usize DEFAULT_DISTANCE_FROM_STR = 1;

  usize mMaxStrUnitLen = DEFAULT_MAX_STR_UNIT_LENGTH;
  usize mMinStrNumUnits = DEFAULT_MIN_STR_UNITS;
  usize mMinStrLength = DEFAULT_MIN_STR_LENGTH;
  usize mDistFromStr = DEFAULT_DISTANCE_FROM_STR;
};

struct StrResult {
  bool mFoundStr = false;
  usize mStrLen = 0;
  std::string mStrMotif;

  friend auto operator==(const StrResult& lhs, const StrResult& rhs) -> bool;
  friend auto operator!=(const StrResult& lhs, const StrResult& rhs) -> bool;
};

[[nodiscard]] auto FindStr(std::string_view seq, usize pos, const StrParams& params = StrParams()) -> StrResult;

#endif  // SRC_LANCET_BASE_FIND_STR_H_
