#ifndef SRC_LANCET_CALLER_RAW_VARIANT_H_
#define SRC_LANCET_CALLER_RAW_VARIANT_H_

#include <string>
#include <utility>

#include "absl/container/flat_hash_map.h"
#include "lancet/base/find_str.h"
#include "lancet/base/types.h"

namespace lancet::caller {

class RawVariant {
 public:
  RawVariant() = default;

  enum class Type : i8 { REF = -1, SNV = 0, INS = 1, DEL = 2, MNP = 3 };
  enum class State : i8 { NONE = -1, SHARED = 0, NORMAL = 1, TUMOR = 2 };

  // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
  usize mChromIndex = -1;
  usize mGenomeStart1 = -1;
  i64 mAlleleLength = -1;
  Type mType = Type::REF;
  std::string mChromName;
  std::string mRefAllele;
  std::string mAltAllele;
  StrResult mStrResult;

  // haplotype index identifier -> start index of variant in haplotype
  absl::flat_hash_map<usize, usize> mHapStart0Idxs;
  // NOLINTEND(misc-non-private-member-variables-in-classes)

  template <typename HashState>
  friend auto AbslHashValue(HashState hash_state, const RawVariant& var) -> HashState {
    return HashState::combine(std::move(hash_state), var.mChromIndex, var.mGenomeStart1, var.mAlleleLength,
                              static_cast<i8>(var.mType), var.mChromName, var.mRefAllele, var.mAltAllele);
  }

  friend auto operator==(const RawVariant& lhs, const RawVariant& rhs) -> bool {
    return lhs.mChromIndex == rhs.mChromIndex && lhs.mGenomeStart1 == rhs.mGenomeStart1 &&
           lhs.mAlleleLength == rhs.mAlleleLength && lhs.mType == rhs.mType && lhs.mChromName == rhs.mChromName &&
           lhs.mRefAllele == rhs.mRefAllele && lhs.mAltAllele == rhs.mAltAllele;
  }

  friend auto operator<(const RawVariant& lhs, const RawVariant& rhs) -> bool {
    // NOLINTBEGIN(readability-braces-around-statements)
    if (lhs.mChromIndex != rhs.mChromIndex) return lhs.mChromIndex < rhs.mChromIndex;
    if (lhs.mGenomeStart1 != rhs.mGenomeStart1) return lhs.mGenomeStart1 < rhs.mGenomeStart1;
    if (lhs.mRefAllele != rhs.mRefAllele) return lhs.mRefAllele < rhs.mRefAllele;
    return lhs.mAltAllele < rhs.mAltAllele;
    // NOLINTEND(readability-braces-around-statements)
  }

  friend auto operator!=(const RawVariant& lhs, const RawVariant& rhs) -> bool { return !(rhs == lhs); }
  friend auto operator>(const RawVariant& lhs, const RawVariant& rhs) -> bool { return rhs < lhs; }
  friend auto operator<=(const RawVariant& lhs, const RawVariant& rhs) -> bool { return !(rhs < lhs); }
  friend auto operator>=(const RawVariant& lhs, const RawVariant& rhs) -> bool { return !(lhs < rhs); }
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_RAW_VARIANT_H_
