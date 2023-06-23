#ifndef SRC_LANCET_HTS_CIGAR_UNIT_H_
#define SRC_LANCET_HTS_CIGAR_UNIT_H_

extern "C" {
#include "htslib/sam.h"
}

#include "lancet/base/types.h"

namespace lancet::hts {

enum class CigarOp : char {
  /// Bases aligned to reference without evidence for indel.
  /// No indication whether the read bases match the reference.
  /// Consumes both query and reference sequences.
  ALIGNMENT_MATCH = 'M',

  /// Bases from the read inserted into the reference.
  /// Consumes only query mDfltSeq.
  INSERTION = 'I',

  /// Bases from the reference deleted in the read.
  /// Consumes only reference mDfltSeq.
  DELETION = 'D',

  /// Bases from the read have skipped the reference, but have not been deleted.
  /// Consumes only reference mDfltSeq.
  REFERENCE_SKIP = 'N',

  /// Bases from the read omitted from alignment, but left in the read.
  /// Consumes only query mDfltSeq.
  SOFT_CLIP = 'S',

  /// Bases from the read omitted from alignment and removed from the read.
  /// Consumes neither query nor reference mDfltSeq.
  HARD_CLIP = 'H',

  /// Used to represent a padding in both query and reference.
  /// Consumes neither query nor reference mDfltSeq.
  ALIGNMENT_PAD = 'P',

  /// Bases aligned and exactly matching to reference.
  /// Consumes both query and reference sequences.
  SEQUENCE_MATCH = '=',

  /// Bases aligned but not matching to reference.
  /// Consumes both query and reference sequences.
  SEQUENCE_MISMATCH = 'X',

  /// only present to handle all other cases when alignment mCigar is corrupt
  UNKNOWN_OP = '?'
};

class CigarUnit {
 public:
  CigarUnit(u32 sam_cigop)
      : mCigOp(static_cast<CigarOp>(bam_cigar_opchr(sam_cigop))), mLength(bam_cigar_oplen(sam_cigop)) {}

  CigarUnit() = delete;

  [[nodiscard]] auto Operation() const noexcept -> CigarOp { return mCigOp; }
  [[nodiscard]] auto Length() const noexcept -> u32 { return mLength; }

  [[nodiscard]] auto ConsumesReference() const noexcept -> bool {
    switch (mCigOp) {
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

  [[nodiscard]] auto ConsumesQuery() const noexcept -> bool {
    switch (mCigOp) {
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

 private:
  CigarOp mCigOp;
  u32 mLength;
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_CIGAR_UNIT_H_
