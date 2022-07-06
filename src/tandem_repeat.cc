#include "lancet2/tandem_repeat.h"

#include <array>
#include "lancet2/assert_macro.h"

namespace lancet2 {
auto FindTandemRepeat(std::string_view seq, std::size_t pos, const TandemRepeatParams& params) -> TandemRepeatResult {
  TandemRepeatResult result{false, 0, ""};

  const auto maxStrUnitLen = static_cast<std::size_t>(params.maxSTRUnitLength);
  const auto minStrUnits = static_cast<std::size_t>(params.minSTRUnitsToReport);
  const auto minStrLen = static_cast<std::size_t>(params.minSTRLengthToReport);
  const auto distFromStr = static_cast<std::size_t>(params.distFromSTR);

  // initialize the offsets
  constexpr std::size_t offsetTableSize = 100;
  std::array<std::array<std::size_t, offsetTableSize>, offsetTableSize> offsets{};

  for (std::size_t merlen = 1; merlen <= maxStrUnitLen; ++merlen) {
    for (std::size_t phase = 0; phase < merlen; ++phase) {
      offsets.at(merlen).at(phase) = phase;
    }
  }

  // now scan the sequence, considering mers starting at position i
  for (std::size_t seqIdx = 0; seqIdx < seq.length(); seqIdx++) {
    // consider all possible merlens from 1 to max
    for (std::size_t merlen = 1; merlen <= maxStrUnitLen; merlen++) {
      const auto phase = seqIdx % merlen;
      const auto offset = offsets.at(merlen).at(phase);

      // compare [i..i+merlen) to [offset..offset+merlen)
      std::size_t endIdx = 0;
      while (endIdx < merlen && seqIdx + endIdx < seq.length() && seq[seqIdx + endIdx] == seq[offset + endIdx]) {
        ++endIdx;
      }

      // is endIdx the end of the tandem?
      if (endIdx != merlen || seqIdx + endIdx + 1 == seq.length()) {
        LANCET_ASSERT(offset + merlen - 1 < seq.length());  // NOLINT

        // am i the leftmost version of this tandem?
        if (offset == 0 || seq[offset - 1] != seq[offset + merlen - 1]) {
          // is it long enough to report?
          if (((seqIdx - offset) / merlen) >= minStrUnits && (seqIdx - offset >= minStrLen)) {
            // is it primitive?
            std::size_t ml = 1;

            while (ml < merlen) {
              const auto units = (seqIdx - offset + endIdx) / ml;

              bool allmatch = true;
              for (std::size_t tmpIdx = 1; allmatch && (tmpIdx < units); tmpIdx++) {
                // compare the bases of the current unit to those of unit0
                for (std::size_t m = 0; m < ml; m++) {
                  if (seq[offset + m] != seq[offset + tmpIdx * ml + m]) {
                    allmatch = false;
                    break;
                  }
                }
              }

              if (!allmatch) {
                ml++;
                continue;
              }

              break;
            }

            // everything checks, now report it
            if (ml == merlen) {
              // start end length
              const auto start = offset > distFromStr ? offset - distFromStr : 0;
              const auto end = seqIdx + endIdx;

              if ((pos >= start) && pos <= (end + distFromStr)) {
                // store STR motif and size
                result.foundSTR = true;
                result.strLength = seqIdx + endIdx - offset;
                for (std::size_t z = 0; z < merlen; ++z) {
                  result.strMotif += seq[offset + z];
                }
              }
            }
          }
        }

        offsets.at(merlen).at(phase) = seqIdx;
      }
    }
  }
  return result;
}
}  // namespace lancet2
