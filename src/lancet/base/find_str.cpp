#include "lancet/base/find_str.h"

#include <array>
#include <string_view>

#include "lancet/base/assert.h"

auto FindStr(std::string_view seq, usize pos, const StrParams& params) -> StrResult {
  StrResult result{.mFoundStr = false, .mStrLen = 0, .mStrMotif = ""};

  // initialize the offsets
  static constexpr usize NUM_OFFSETS = 100;
  std::array<std::array<usize, NUM_OFFSETS>, NUM_OFFSETS> offsets{};

  for (usize merlen = 1; merlen <= params.mMaxStrUnitLen; ++merlen) {
    for (usize phase = 0; phase < merlen; ++phase) {
      offsets.at(merlen).at(phase) = phase;
    }
  }

  // now scan the sequence, considering mers starting at position idx
  for (usize bpos = 0; bpos < seq.length(); bpos++) {
    // consider all possible merlens from 1 to max_str_unit_len
    for (usize merlen = 1; merlen <= params.mMaxStrUnitLen; merlen++) {
      const auto phase = bpos % merlen;
      const auto offset = offsets.at(merlen).at(phase);

      // compare [i..i+merlen) to [offset..offset+merlen)
      usize end_idx = 0;
      while (end_idx < merlen && bpos + end_idx < seq.length() && seq[bpos + end_idx] == seq[offset + end_idx]) {
        ++end_idx;
      }

      // is end_idx the end of the tandem ?
      if (end_idx != merlen || bpos + end_idx + 1 == seq.length()) {
        LANCET_ASSERT(offset + merlen - 1 < seq.length())
        // am i the leftmost version of this tandem ?
        if (offset == 0 || seq[offset - 1] != seq[offset + merlen - 1]) {
          // is it long enough to report ?
          if (((bpos - offset) / merlen) >= params.mMinStrNumUnits && (bpos - offset >= params.mMinStrLength)) {
            // is it primitive?
            usize mlen = 1;
            while (mlen < merlen) {
              const auto units = (bpos - offset + end_idx) / mlen;

              bool allmatch = true;
              for (usize tmp_idx = 1; allmatch && (tmp_idx < units); tmp_idx++) {
                // compare the bases of the current unit to those of unit0
                for (usize other = 0; other < mlen; other++) {
                  if (seq[offset + other] != seq[offset + tmp_idx * mlen + other]) {
                    allmatch = false;
                    break;
                  }
                }
              }

              if (!allmatch) {
                mlen++;
                continue;
              }

              break;
            }

            // everything checks, now report it
            if (mlen == merlen) {
              // start end length
              const auto start = offset > params.mDistFromStr ? offset - params.mDistFromStr : 0;
              const auto end = bpos + end_idx;

              if ((pos >= start) && pos <= (end + params.mDistFromStr)) {
                // store STR motif and size
                result.mFoundStr = true;
                result.mStrLen = bpos + end_idx - offset;
                for (usize midx = 0; midx < merlen; ++midx) {
                  result.mStrMotif += seq[offset + midx];
                }
              }
            }
          }
        }

        offsets.at(merlen).at(phase) = bpos;
      }
    }
  }

  return result;
}

auto operator==(const StrResult& lhs, const StrResult& rhs) -> bool {
  return lhs.mFoundStr == rhs.mFoundStr && lhs.mStrLen == rhs.mStrLen && lhs.mStrMotif == rhs.mStrMotif;
}

auto operator!=(const StrResult& lhs, const StrResult& rhs) -> bool { return !(rhs == lhs); }
