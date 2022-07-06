#include "lancet2/cigar.h"

#include "absl/strings/str_format.h"

namespace lancet2 {
auto CigarUnit::ConsumesReference() const -> bool {
  switch (Operation) {
    case CigarOp::ALIGNMENT_MATCH:
    case CigarOp::DELETION:
    case CigarOp::REFERENCE_SKIP:
    case CigarOp::SEQUENCE_MATCH:
    case CigarOp::SEQUENCE_MISMATCH:
      return true;

    case CigarOp::INSERTION:
    case CigarOp::SOFT_CLIP:
    case CigarOp::HARD_CLIP:
    case CigarOp::ALIGNMENT_PAD:
    case CigarOp::UNKNOWN_OP:
    default:
      return false;
  }
}

auto CigarUnit::ConsumesQuery() const -> bool {
  switch (Operation) {
    case CigarOp::ALIGNMENT_MATCH:
    case CigarOp::INSERTION:
    case CigarOp::SOFT_CLIP:
    case CigarOp::SEQUENCE_MATCH:
    case CigarOp::SEQUENCE_MISMATCH:
      return true;

    case CigarOp::DELETION:
    case CigarOp::REFERENCE_SKIP:
    case CigarOp::HARD_CLIP:
    case CigarOp::ALIGNMENT_PAD:
    case CigarOp::UNKNOWN_OP:
    default:
      return false;
  }
}

auto CigarUnit::ToString() const -> std::string {
  return absl::StrFormat("%d%c", Length, static_cast<char>(Operation));
}
}  // namespace lancet2
