#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include "lancet2/core_enums.h"

namespace lancet2 {
struct ReadInfo {
  std::string readName;                     // NOLINT
  std::string chromName;                    // NOLINT
  std::string sequence;                     // NOLINT
  std::string quality;                      // NOLINT
  std::int64_t startPos0 = -1;              // NOLINT
  Strand strand = Strand::FWD;              // NOLINT
  SampleLabel label = SampleLabel::NORMAL;  // NOLINT
  std::int8_t haplotypeID = -1;             // NOLINT
  std::string tenxBarcode;                  // NOLINT

  [[nodiscard]] auto Length() const noexcept -> std::size_t { return sequence.length(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool {
    return sequence.empty() && quality.empty() && chromName.empty() && startPos0 == -1;
  }
};
}  // namespace lancet2
