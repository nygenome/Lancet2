#ifndef SRC_LANCET_HTS_CIGAR_UTILS_H_
#define SRC_LANCET_HTS_CIGAR_UTILS_H_

#include <vector>

#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/hts/cigar_unit.h"

namespace lancet::hts {

// ============================================================================
// ComputeEditDistance: SAM specification NM tag computation.
//
// NM = edit distance to the reference = mismatches + insertions + deletions.
// Since minimap2 emits the M (alignment match) CIGAR operator rather than
// the newer =/X (sequence match/mismatch) operators, mismatches under M
// blocks cannot be determined from the CIGAR alone — the actual encoded
// query and target sequences are compared base-by-base.
//
// Per SAM spec (https://samtools.github.io/hts-specs/SAMtags.pdf):
//   M  — compare query[qpos] vs target[tpos], +1 if different
//   =  — sequence match, no NM contribution
//   X  — sequence mismatch, +len to NM
//   I  — insertion, +len to NM
//   D  — deletion, +len to NM
//   S  — soft clip, advance query only, excluded from NM
//   N  — reference skip, advance reference only, excluded from NM
//   H,P — no advancement, excluded from NM
// ============================================================================
[[nodiscard]] inline auto ComputeEditDistance(const std::vector<CigarUnit>& cigar,
                                              absl::Span<const u8> encoded_query,
                                              absl::Span<const u8> encoded_target) -> u32 {
  u32 nm = 0;
  usize qpos = 0;
  usize tpos = 0;

  for (const auto& unit : cigar) {
    const auto op = unit.Operation();
    const u32 len = unit.Length();

    switch (op) {
      case CigarOp::ALIGNMENT_MATCH:
        // M op: must compare sequences base-by-base for mismatches.
        // Encoded bases are numeric (0-4), so comparison is case-insensitive.
        for (u32 i = 0; i < len; ++i, ++qpos, ++tpos) {
          if (qpos < encoded_query.size() && tpos < encoded_target.size()
              && encoded_query[qpos] != encoded_target[tpos]) {
            ++nm;
          }
        }
        break;
      case CigarOp::SEQUENCE_MATCH:
        // = op: all matches, no NM contribution
        qpos += len;
        tpos += len;
        break;
      case CigarOp::SEQUENCE_MISMATCH:
        // X op: all mismatches
        nm += len;
        qpos += len;
        tpos += len;
        break;
      case CigarOp::INSERTION:
        nm += len;
        qpos += len;
        break;
      case CigarOp::DELETION:
        nm += len;
        tpos += len;
        break;
      case CigarOp::SOFT_CLIP:
        qpos += len;
        break;
      case CigarOp::REFERENCE_SKIP:
        tpos += len;
        break;
      default:
        break;  // H, P: no advancement, no NM
    }
  }

  return nm;
}

// ============================================================================
// CigarRefPosToQueryPos: map a reference position to a query position.
//
// Walks the CIGAR consuming reference and query positions as appropriate.
// Returns the query offset at the first alignment column covering `ref_pos`.
// If ref_pos falls in a deletion (no query base), returns the query position
// at the start of that deletion.
// ============================================================================
[[nodiscard]] inline auto CigarRefPosToQueryPos(const std::vector<CigarUnit>& cigar,
                                                 const usize ref_pos) -> usize {
  usize qpos = 0;
  usize tpos = 0;

  for (const auto& unit : cigar) {
    const auto op = unit.Operation();
    const u32 len = unit.Length();

    switch (op) {
      case CigarOp::ALIGNMENT_MATCH:
      case CigarOp::SEQUENCE_MATCH:
      case CigarOp::SEQUENCE_MISMATCH:
        for (u32 i = 0; i < len; ++i, ++qpos, ++tpos) {
          // NOLINTNEXTLINE(readability-braces-around-statements)
          if (tpos == ref_pos) return qpos;
        }
        break;
      case CigarOp::INSERTION:
        qpos += len;
        break;
      case CigarOp::DELETION:
      case CigarOp::REFERENCE_SKIP:
        for (u32 i = 0; i < len; ++i, ++tpos) {
          // NOLINTNEXTLINE(readability-braces-around-statements)
          if (tpos == ref_pos) return qpos;
        }
        break;
      case CigarOp::SOFT_CLIP:
        qpos += len;
        break;
      default:
        break;  // H, P
    }
  }

  return qpos;  // ref_pos beyond CIGAR — return end of query
}

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_CIGAR_UTILS_H_
