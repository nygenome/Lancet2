#include "lancet2/hts_alignment.h"

#include "absl/strings/ascii.h"
#include "lancet2/assert_macro.h"

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-align"
#pragma clang diagnostic ignored "-Wcast-qual"
#elif defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif

#include "htslib/sam.h"

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace lancet2 {
auto HtsAlignment::BuildReadInfo(SampleLabel label, u8 min_bq, u8 max_kmer_size) const -> ReadInfo {
  const auto seqLen = readSequence.length();
  const auto qualLen = readQuality.length();
  LANCET_ASSERT(seqLen == qualLen);  // NOLINT

  usize trim5 = 0;
  for (trim5 = 0; trim5 < seqLen; ++trim5) {
    const auto base = absl::ascii_toupper(static_cast<unsigned char>(readSequence[trim5]));
    if ((base == 'A' || base == 'C' || base == 'G' || base == 'T') && static_cast<int>(readQuality[trim5]) >= min_bq) {
      break;
    }
  }

  // return empty read info
  if (trim5 == seqLen) return ReadInfo{};

  usize trim3 = 0;
  for (auto idx = seqLen - 1; idx >= 0; --idx) {
    const auto base = absl::ascii_toupper(static_cast<unsigned char>(readSequence[idx]));
    if ((base == 'A' || base == 'C' || base == 'G' || base == 'T') && static_cast<int>(readQuality[idx]) >= min_bq) {
      break;
    }
    trim3++;
  }

  if ((seqLen - trim5 - trim3) < static_cast<usize>(max_kmer_size)) return ReadInfo{};

  ReadInfo ri;
  ri.readName = readName;
  ri.chromName = contig;
  ri.sequence = readSequence.substr(trim5, seqLen - trim5 - trim3);
  ri.quality = readQuality.substr(trim5, qualLen - trim5 - trim3);
  ri.startPos0 = startPosition0 + static_cast<i64>(trim5);
  ri.strand = GetReadStrand();
  ri.label = label;

  const auto hpItr = tagsData.find("HP");
  if (hpItr != tagsData.end()) {
    ri.haplotypeID = static_cast<i8>(bam_aux2i(hpItr->second));
  }

  const auto bxItr = tagsData.find("BX");
  if (bxItr != tagsData.end()) {
    ri.tenxBarcode = bam_aux2Z(bxItr->second);
  }

  return ri;
}

auto HtsAlignment::GetReadStrand() const -> Strand { return IsReverseStrand() ? Strand::REV : Strand::FWD; }
auto HtsAlignment::IsDuplicate() const -> bool { return (samFlags & BAM_FDUP) != 0; }                // NOLINT
auto HtsAlignment::IsSupplementary() const -> bool { return (samFlags & BAM_FSUPPLEMENTARY) != 0; }  // NOLINT
auto HtsAlignment::IsPrimary() const -> bool { return (samFlags & BAM_FSECONDARY) == 0; }            // NOLINT
auto HtsAlignment::IsSecondary() const -> bool { return (samFlags & BAM_FSECONDARY) != 0; }          // NOLINT
auto HtsAlignment::IsQcFailed() const -> bool { return (samFlags & BAM_FQCFAIL) != 0; }              // NOLINT
auto HtsAlignment::IsUnmapped() const -> bool { return (samFlags & BAM_FUNMAP) != 0; }               // NOLINT
auto HtsAlignment::IsMateUnmapped() const -> bool { return (samFlags & BAM_FMUNMAP) != 0; }          // NOLINT
auto HtsAlignment::IsReverseStrand() const -> bool { return (samFlags & BAM_FREVERSE) != 0; }        // NOLINT
auto HtsAlignment::IsMateReverseStrand() const -> bool { return (samFlags & BAM_FMREVERSE) != 0; }   // NOLINT
auto HtsAlignment::IsPaired() const -> bool { return (samFlags & BAM_FPAIRED) != 0; }                // NOLINT
auto HtsAlignment::IsProperPair() const -> bool { return (samFlags & BAM_FPROPER_PAIR) != 0; }       // NOLINT
auto HtsAlignment::IsRead1() const -> bool { return (samFlags & BAM_FREAD1) != 0; }                  // NOLINT
auto HtsAlignment::IsRead2() const -> bool { return (samFlags & BAM_FREAD2) != 0; }                  // NOLINT

auto HtsAlignment::GetSoftClips(std::vector<u32> *clip_sizes, std::vector<u32> *read_positions,
                                std::vector<u32> *genome_positions, bool use_padded) const -> bool {
  // initialize positions & flags
  auto refPosition = static_cast<u32>(startPosition0);
  u32 readPosition = 0;
  bool softClipFound = false;
  bool firstCigarOp = true;

  for (const auto &cigUnit : cigar) {
    switch (cigUnit.Operation) {
      // increase both read & genome positions on CIGAR chars [DMXN=]
      case CigarOp::DELETION:
      case CigarOp::ALIGNMENT_MATCH:
      case CigarOp::SEQUENCE_MISMATCH:
      case CigarOp::REFERENCE_SKIP:
      case CigarOp::SEQUENCE_MATCH:
        refPosition += cigUnit.Length;
        readPosition += cigUnit.Length;
        break;

        // increase read position on insertion, genome position only if @usePadded is true
      case CigarOp::INSERTION:
        readPosition += cigUnit.Length;
        if (use_padded) refPosition += cigUnit.Length;
        break;

      case CigarOp::SOFT_CLIP:
        softClipFound = true;

        //////////////////////////////////////////////////////////////////////////////
        // if we are dealing with the *first* CIGAR operation
        // for this alignment, we increment the read position so that
        // the read and genome position of the clip are referring to the same base.
        // For example, in the alignment below, the ref position would be 4, yet
        //              the read position would be 0. Thus, to "sync" the two,
        //              we need to increment the read position by the length of the
        //              soft clip.
        // Read:  ATCGTTTCGTCCCTGC
        // Ref:   GGGATTTCGTCCCTGC
        // Cigar: SSSSMMMMMMMMMMMM
        //
        // NOTE: This only needs to be done if the soft clip is the _first_ CIGAR op.
        //////////////////////////////////////////////////////////////////////////////
        if (firstCigarOp) readPosition += cigUnit.Length;

        // track the soft clip's size, read position, and genome position
        if (clip_sizes != nullptr) clip_sizes->push_back(cigUnit.Length);
        if (read_positions != nullptr) read_positions->push_back(readPosition);
        if (genome_positions != nullptr) genome_positions->push_back(refPosition);
        break;

        // any other CIGAR operations have no effect
      case CigarOp::HARD_CLIP:
      case CigarOp::ALIGNMENT_PAD:
      case CigarOp::UNKNOWN_OP:
      default:
        break;
    }

    // clear our "first pass" flag
    firstCigarOp = false;
  }

  return softClipFound;
}
}  // namespace lancet2
