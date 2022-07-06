#pragma once

#include <string>

#include "lancet2/core_enums.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
struct ReadInfo {
  std::string readName;                     // NOLINT
  std::string chromName;                    // NOLINT
  std::string sequence;                     // NOLINT
  std::string quality;                      // NOLINT
  i64 startPos0 = -1;                       // NOLINT
  Strand strand = Strand::FWD;              // NOLINT
  SampleLabel label = SampleLabel::NORMAL;  // NOLINT
  i8 haplotypeID = -1;                      // NOLINT
  std::string tenxBarcode;                  // NOLINT

  [[nodiscard]] auto Length() const noexcept -> usize { return sequence.length(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool {
    return sequence.empty() && quality.empty() && chromName.empty() && startPos0 == -1;
  }
};
}  // namespace lancet2
