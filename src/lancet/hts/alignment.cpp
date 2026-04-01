#include "lancet/hts/alignment.h"

#include <algorithm>
#include <array>
#include <string>
#include <string_view>
#include <vector>

#include "lancet/hts/cigar_unit.h"

extern "C" {
#include "htslib/sam.h"
}

#include "absl/strings/str_cat.h"
#include "lancet/base/types.h"
#include "lancet/hts/reference.h"
#include "spdlog/fmt/bundled/core.h"

namespace lancet::hts {

// 4-bit BAM nucleotide encoding -> ASCII character lookup table
static constexpr auto SEQ_4BIT_TO_CHAR =
    std::array<char, 16>{'N', 'A', 'C', 'N', 'G', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};

// ---------------------------------------------------------------------------
// Lifecycle management
// ---------------------------------------------------------------------------

void Alignment::ClearAllFields() {
  mStart0 = -1;
  mMateStart0 = -1;
  mInsertSize = -1;
  mChromIdx = -1;
  mMateChromIdx = -1;
  mSamFlag = 0;
  mMapQual = 0;
  mRawAln = nullptr;
}

void Alignment::PopulateFromRaw(bam1_t* aln) {
  // Cache the cheap scalar core fields to avoid repeated pointer chasing
  mStart0 = aln->core.pos;
  mMateStart0 = aln->core.mpos;
  mInsertSize = aln->core.isize;
  mChromIdx = aln->core.tid;
  mMateChromIdx = aln->core.mtid;
  mSamFlag = aln->core.flag;
  mMapQual = aln->core.qual;

  // Retain the raw pointer for zero-copy field access.
  // WARNING: This pointer is managed by the Iterator and becomes invalid on ++itr.
  mRawAln = aln;
}

// ---------------------------------------------------------------------------
// Zero-copy proxy methods
// ---------------------------------------------------------------------------

auto Alignment::QnameView() const noexcept -> std::string_view {
  if (mRawAln == nullptr) return {};
  return bam_get_qname(mRawAln);
}

auto Alignment::Length() const noexcept -> usize {
  if (mRawAln == nullptr) return 0;
  return static_cast<usize>(mRawAln->core.l_qseq);
}

auto Alignment::CigarData() const -> std::vector<CigarUnit> {
  if (mRawAln == nullptr) return {};
  const auto* raw_cigar = bam_get_cigar(mRawAln);
  const auto n_cigar = static_cast<usize>(mRawAln->core.n_cigar);
  std::vector<CigarUnit> result;
  result.reserve(n_cigar);
  for (usize idx = 0; idx < n_cigar; ++idx) {
    result.emplace_back(CigarUnit(raw_cigar[idx]));
  }
  return result;
}

auto Alignment::CigarString() const -> std::string {
  if (mRawAln == nullptr) return {};
  const auto* raw_cigar = bam_get_cigar(mRawAln);
  const auto n_cigar = static_cast<usize>(mRawAln->core.n_cigar);
  std::string result;
  result.reserve(n_cigar * 4);
  for (usize idx = 0; idx < n_cigar; ++idx) {
    absl::StrAppend(&result, bam_cigar_oplen(raw_cigar[idx]), std::string(1, bam_cigar_opchr(raw_cigar[idx])));
  }
  result.shrink_to_fit();
  return result;
}

// ---------------------------------------------------------------------------
// On-demand deep extraction methods
// ---------------------------------------------------------------------------

auto Alignment::BuildSequence() const -> std::string {
  if (mRawAln == nullptr) return {};
  const auto seq_len = static_cast<usize>(mRawAln->core.l_qseq);
  const auto* raw_seq = bam_get_seq(mRawAln);
  std::string result(seq_len, 'N');
  for (usize idx = 0; idx < seq_len; ++idx) {
    result[idx] = SEQ_4BIT_TO_CHAR[bam_seqi(raw_seq, idx)];
  }
  return result;
}

auto Alignment::BuildQualities() const -> std::vector<u8> {
  if (mRawAln == nullptr) return {};
  const auto seq_len = static_cast<usize>(mRawAln->core.l_qseq);
  const auto* raw_qual = bam_get_qual(mRawAln);
  return {raw_qual, raw_qual + seq_len};
}

// ---------------------------------------------------------------------------
// Aux tag access: direct bam_aux_get routing
// ---------------------------------------------------------------------------

auto Alignment::FindRawTag(std::string_view tag_name) const noexcept -> const u8* {
  if (mRawAln == nullptr || tag_name.size() != 2) return nullptr;
  // bam_aux_get expects a two-character null-terminated string
  const char tag_str[3] = {tag_name[0], tag_name[1], '\0'};
  return bam_aux_get(mRawAln, tag_str);
}

auto Alignment::HasTag(std::string_view tag_name) const noexcept -> bool { return FindRawTag(tag_name) != nullptr; }

// Explicit template specializations for ExtractTagValue
template <>
auto Alignment::ExtractTagValue<i64>(const u8* raw_aux, std::string_view tag_name) -> absl::StatusOr<i64> {
  const auto val = bam_aux2i(raw_aux);
  if (errno == EINVAL) {
    return absl::Status(absl::StatusCode::kInvalidArgument,
                        fmt::format("Tag {} does not contain an integer value", tag_name));
  }
  return val;
}

template <>
auto Alignment::ExtractTagValue<f64>(const u8* raw_aux, std::string_view tag_name) -> absl::StatusOr<f64> {
  const auto val = bam_aux2f(raw_aux);
  if (errno == EINVAL) {
    return absl::Status(absl::StatusCode::kInvalidArgument,
                        fmt::format("Tag {} does not contain a float value", tag_name));
  }
  return val;
}

template <>
auto Alignment::ExtractTagValue<std::string_view>(const u8* raw_aux, std::string_view tag_name)
    -> absl::StatusOr<std::string_view> {
  const auto* val = bam_aux2Z(raw_aux);
  if (val == nullptr) {
    return absl::Status(absl::StatusCode::kInvalidArgument,
                        fmt::format("Tag {} does not contain a string value", tag_name));
  }
  return std::string_view(val);
}

// ---------------------------------------------------------------------------
// BitwiseFlag implementation
// ---------------------------------------------------------------------------

auto Alignment::BitwiseFlag::GetStrand() const noexcept -> Strand { return IsFwdStrand() ? Strand::FWD : Strand::REV; }
auto Alignment::BitwiseFlag::GetMateStrand() const noexcept -> Strand {
  return IsMateFwdStrand() ? Strand::FWD : Strand::REV;
}
auto Alignment::BitwiseFlag::IsFwdStrand() const noexcept -> bool { return (mFlag & BAM_FREVERSE) == 0; }
auto Alignment::BitwiseFlag::IsRevStrand() const noexcept -> bool { return (mFlag & BAM_FREVERSE) != 0; }
auto Alignment::BitwiseFlag::IsMateFwdStrand() const noexcept -> bool { return (mFlag & BAM_FMREVERSE) == 0; }
auto Alignment::BitwiseFlag::IsMateRevStrand() const noexcept -> bool { return (mFlag & BAM_FMREVERSE) != 0; }
auto Alignment::BitwiseFlag::IsQcFail() const noexcept -> bool { return (mFlag & BAM_FQCFAIL) != 0; }
auto Alignment::BitwiseFlag::IsDuplicate() const noexcept -> bool { return (mFlag & BAM_FDUP) != 0; }
auto Alignment::BitwiseFlag::IsPrimary() const noexcept -> bool { return (mFlag & BAM_FSECONDARY) == 0; }
auto Alignment::BitwiseFlag::IsSecondary() const noexcept -> bool { return (mFlag & BAM_FSECONDARY) != 0; }
auto Alignment::BitwiseFlag::IsSupplementary() const noexcept -> bool { return (mFlag & BAM_FSUPPLEMENTARY) != 0; }
auto Alignment::BitwiseFlag::IsMapped() const noexcept -> bool { return (mFlag & BAM_FUNMAP) == 0; }
auto Alignment::BitwiseFlag::IsUnmapped() const noexcept -> bool { return (mFlag & BAM_FUNMAP) != 0; }
auto Alignment::BitwiseFlag::IsMateMapped() const noexcept -> bool { return (mFlag & BAM_FMUNMAP) == 0; }
auto Alignment::BitwiseFlag::IsMateUnmapped() const noexcept -> bool { return (mFlag & BAM_FMUNMAP) != 0; }
auto Alignment::BitwiseFlag::IsPairedInSequencing() const noexcept -> bool { return (mFlag & BAM_FPAIRED) != 0; }
auto Alignment::BitwiseFlag::IsMappedProperPair() const noexcept -> bool { return (mFlag & BAM_FPROPER_PAIR) != 0; }
auto Alignment::BitwiseFlag::IsRead1() const noexcept -> bool { return (mFlag & BAM_FREAD1) != 0; }
auto Alignment::BitwiseFlag::IsRead2() const noexcept -> bool { return (mFlag & BAM_FREAD2) != 0; }
auto Alignment::BitwiseFlag::HasFlagsSet(u16 check_flags) const noexcept -> bool { return (mFlag & check_flags) != 0; }
auto Alignment::BitwiseFlag::HasFlagsUnset(u16 check_flags) const noexcept -> bool {
  return (mFlag & check_flags) == 0;
}

// ---------------------------------------------------------------------------
// Location and region helpers
// ---------------------------------------------------------------------------

auto Alignment::MateLocation() const noexcept -> MateInfo {
  return {.mChromIndex = mMateChromIdx, .mMateStartPos0 = mMateStart0};
}

auto Alignment::MateOverlapsRegion(const Reference::Region& region) const noexcept -> bool {
  return region.ChromIndex() == mMateChromIdx && (mMateStart0 + 1) >= region.StartPos1() &&
         (mMateStart0 + 1) <= region.EndPos1();
}

auto Alignment::OverlapsRegion(const Reference::Region& region) const noexcept -> bool {
  return mChromIdx == region.ChromIndex() && (mStart0 + 1) >= region.StartPos1() && (mStart0 + 1) <= region.EndPos1();
}

auto Alignment::IsEmpty() const noexcept -> bool {
  return mStart0 == -1 && mMateStart0 == -1 && mInsertSize == -1 && mChromIdx == -1 && mMateChromIdx == -1 &&
         mSamFlag == 0 && mMapQual == 0 && mRawAln == nullptr;
}

// ---------------------------------------------------------------------------
// ToString (debug/diagnostic output)
// ---------------------------------------------------------------------------

auto Alignment::ToString(const Reference& ref) const -> std::string {
  const auto chrom = ref.FindChromByIndex(mChromIdx);
  const auto mate_chrom = ref.FindChromByIndex(mMateChromIdx);

  const auto both_chroms_same = chrom.ok() && mate_chrom.ok() && chrom->Name() == mate_chrom->Name();
  // NOLINTNEXTLINE(readability-avoid-nested-conditional-operator)
  const auto rnext = both_chroms_same ? "=" : mate_chrom.ok() ? mate_chrom->Name() : "*";

  const auto seq = BuildSequence();
  const auto quals = BuildQualities();

  std::string fastq_quality;
  fastq_quality.resize(quals.size());
  static constexpr u8 PHRED_QUALITY_OFFSET = 33;
  std::transform(quals.cbegin(), quals.cend(), fastq_quality.begin(),
                 [](const u8 base_qual) { return static_cast<char>(base_qual + PHRED_QUALITY_OFFSET); });

  const auto qname_sv = QnameView();

  return fmt::format(
      "{QNAME}\t{FLAG}\t{RNAME}\t{POS}\t{MAPQ}\t{CIGAR}\t{RNEXT}\t{PNEXT}\t{TLEN}\t{SEQ}\t{QUAL}\n",
      fmt::arg("QNAME", qname_sv.empty() ? "*" : std::string(qname_sv)), fmt::arg("FLAG", mSamFlag),
      fmt::arg("RNAME", chrom.ok() ? chrom->Name() : "*"), fmt::arg("POS", mStart0 >= 0 ? mStart0 + 1 : 0),
      fmt::arg("MAPQ", mMapQual), fmt::arg("CIGAR", CigarString()), fmt::arg("RNEXT", rnext),
      fmt::arg("PNEXT", mMateStart0 >= 0 ? mMateStart0 + 1 : 0), fmt::arg("TLEN", mInsertSize), fmt::arg("SEQ", seq),
      fmt::arg("QUAL", fastq_quality));
}

// ---------------------------------------------------------------------------
// Equality operators
// ---------------------------------------------------------------------------

auto Alignment::operator==(const Alignment& rhs) const -> bool {
  return mStart0 == rhs.mStart0 && mMateStart0 == rhs.mMateStart0 && mInsertSize == rhs.mInsertSize &&
         mChromIdx == rhs.mChromIdx && mMateChromIdx == rhs.mMateChromIdx && mSamFlag == rhs.mSamFlag &&
         mMapQual == rhs.mMapQual;
}

auto Alignment::operator!=(const Alignment& rhs) const -> bool { return !(rhs == *this); }

// ---------------------------------------------------------------------------
// Soft-clip detection (reads directly from bam1_t cigar)
// ---------------------------------------------------------------------------

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto Alignment::GetSoftClips(std::vector<u32>* clip_sizes, std::vector<u32>* read_positions,
                             std::vector<u32>* genome_positions, bool use_padded) const -> bool {
  if (mRawAln == nullptr) return false;

  const auto* raw_cigar = bam_get_cigar(mRawAln);
  const auto n_cigar = static_cast<usize>(mRawAln->core.n_cigar);

  // initialize positions & flags
  auto ref_position = static_cast<u32>(mStart0);
  u32 read_position = 0;
  bool soft_clip_found = false;
  bool first_cigar_op = true;

  for (usize idx = 0; idx < n_cigar; ++idx) {
    const auto cig_unit = CigarUnit(raw_cigar[idx]);
    switch (cig_unit.Operation()) {
      // increase both read & genome positions on CIGAR chars [DMXN=]
      case CigarOp::DELETION:
      case CigarOp::ALIGNMENT_MATCH:
      case CigarOp::SEQUENCE_MISMATCH:
      case CigarOp::REFERENCE_SKIP:
      case CigarOp::SEQUENCE_MATCH:
        ref_position += cig_unit.Length();
        read_position += cig_unit.Length();
        break;

        // increase read position on insertion, genome position only if @usePadded is true
      case CigarOp::INSERTION:
        read_position += cig_unit.Length();
        // NOLINTNEXTLINE(readability-braces-around-statements)
        if (use_padded) ref_position += cig_unit.Length();
        break;

      case CigarOp::SOFT_CLIP:
        soft_clip_found = true;

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
        // NOLINTBEGIN(readability-braces-around-statements)
        if (first_cigar_op) read_position += cig_unit.Length();

        // track the soft clip's size, read position, and genome position
        if (clip_sizes != nullptr) clip_sizes->push_back(cig_unit.Length());
        if (read_positions != nullptr) read_positions->push_back(read_position);
        if (genome_positions != nullptr) genome_positions->push_back(ref_position);
        // NOLINTEND(readability-braces-around-statements)
        break;

        // any other CIGAR operations have no effect
      case CigarOp::HARD_CLIP:
      case CigarOp::ALIGNMENT_PAD:
      case CigarOp::UNKNOWN_OP:
      default:
        break;
    }

    // clear our "first pass" flag
    first_cigar_op = false;
  }

  return soft_clip_found;
}

}  // namespace lancet::hts
