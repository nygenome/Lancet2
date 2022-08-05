#pragma once

#include <vector>

#include "absl/types/span.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {

struct BaseCov {
  u32 FwdRaw = 0;
  u32 RevRaw = 0;
  u32 FwdBQPass = 0;
  u32 RevBQPass = 0;
};

using KmerCov = std::vector<BaseCov>;

}  // namespace lancet2
