#pragma once

#include <string>
#include <vector>

#include "lancet2/sized_ints.h"

namespace lancet2 {
enum class CigarOp : char {
  /// Bases aligned to reference without evidence for indel.
  /// No indication whether the read bases match the reference.
  /// Consumes both query and reference sequences.
  ALIGNMENT_MATCH = 'M',

  /// Bases from the read inserted into the reference.
  /// Consumes only query sequence.
  INSERTION = 'I',

  /// Bases from the reference deleted in the read.
  /// Consumes only reference sequence.
  DELETION = 'D',

  /// Bases from the read have skipped the reference, but have not been deleted.
  /// Consumes only reference sequence.
  REFERENCE_SKIP = 'N',

  /// Bases from the read omitted from alignment, but left in the read.
  /// Consumes only query sequence.
  SOFT_CLIP = 'S',

  /// Bases from the read omitted from alignment and removed from the read.
  /// Consumes neither query nor reference sequence.
  HARD_CLIP = 'H',

  /// Used to represent a padding in both query and reference.
  /// Consumes neither query nor reference sequence.
  ALIGNMENT_PAD = 'P',

  /// Bases aligned and exactly matching to reference.
  /// Consumes both query and reference sequences.
  SEQUENCE_MATCH = '=',

  /// Bases aligned but not matching to reference.
  /// Consumes both query and reference sequences.
  SEQUENCE_MISMATCH = 'X',

  /// only present to handle all other cases when alignment cigar is corrupt
  UNKNOWN_OP = '?'
};

class CigarUnit {
 public:
  CigarOp Operation;  // NOLINT
  u32 Length;         // NOLINT

  CigarUnit() = delete;
  explicit CigarUnit(CigarOp op, u32 len) : Operation(op), Length(len) {}
  explicit CigarUnit(const char& op_char, u32 len) {
    switch (op_char) {
      case 'M':
      case 'I':
      case 'D':
      case 'N':
      case 'S':
      case 'H':
      case 'P':
      case '=':
      case 'X':
        Operation = static_cast<CigarOp>(op_char);
        Length = len;
        break;

      default:
        Operation = CigarOp::UNKNOWN_OP;
        Length = len;
    }
  }

  [[nodiscard]] auto ConsumesReference() const -> bool;
  [[nodiscard]] auto ConsumesQuery() const -> bool;

  [[nodiscard]] auto ToString() const -> std::string;
};

using AlignmentCigar = std::vector<CigarUnit>;
}  // namespace lancet2
