#include "lancet/hts/alignment.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <stdexcept>

#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "spdlog/fmt/fmt.h"

namespace lancet::hts {

static constexpr auto SEQ_4BIT_TO_CHAR =
    std::array<char, 16>{'N', 'A', 'C', 'N', 'G', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};

namespace {

inline auto SequenceFrom4Bit(absl::Span<const u8> bases) -> std::string {
  std::string result(bases.length(), 'N');
  for (usize idx = 0; idx < bases.length(); ++idx) {
    result[idx] = SEQ_4BIT_TO_CHAR.at(bam_seqi(bases.data(), idx));
  }
  return result;
}

}  // namespace

void Alignment::ClearAllFields() {
  mStart0 = -1;
  mMateStart0 = -1;
  mInsertSize = -1;
  mChromIdx = -1;
  mMateChromIdx = -1;
  mSamFlag = 0;
  mMapQual = 0;
  mQname.clear();
  mSeq.clear();
  mQual.clear();
  mCigar.clear();
  mAuxTags.clear();
}

void Alignment::PopulateRequestedFields(bam1_t* aln, const Fields fields, const TagNamesSet* fill_tags) {
  PopulateCoreQname(aln);
  if (fields == Alignment::Fields::CORE_QNAME) return;  // NOLINT(readability-braces-around-statements)

  PopulateCigarSeqQual(aln);
  if (fields == Alignment::Fields::CIGAR_SEQ_QUAL) return;  // NOLINT(readability-braces-around-statements)

  return PopulateAuxRgAux(aln, fill_tags);
}

void Alignment::PopulateCoreQname(bam1_t* aln) {
  mStart0 = aln->core.pos;
  mMateStart0 = aln->core.mpos;
  mInsertSize = aln->core.isize;
  mChromIdx = aln->core.tid;
  mMateChromIdx = aln->core.mtid;
  mSamFlag = aln->core.flag;
  mMapQual = aln->core.qual;
  mQname.assign(bam_get_qname(aln));
}

void Alignment::PopulateCigarSeqQual(bam1_t* aln) {
  const absl::Span<const u8> raw_bases = absl::MakeConstSpan(bam_get_seq(aln), aln->core.l_qseq);
  mSeq = SequenceFrom4Bit(raw_bases);

  const absl::Span<const u8> raw_quals = absl::MakeConstSpan(bam_get_qual(aln), aln->core.l_qseq);
  mQual = std::vector<u8>(raw_quals.cbegin(), raw_quals.cend());

  const absl::Span<const u32> raw_cigar = absl::MakeConstSpan(bam_get_cigar(aln), aln->core.n_cigar);
  mCigar = std::vector<u32>(raw_cigar.cbegin(), raw_cigar.cend());
}

void Alignment::PopulateAuxRgAux(bam1_t* aln, const TagNamesSet* fill_tags) {
  if (fill_tags->empty()) return;  // NOLINT(readability-braces-around-statements)

  mAuxTags.clear();
  mAuxTags.reserve(fill_tags->size());

  u8* curr_aux = bam_aux_first(aln);
  while (curr_aux != nullptr) {
    if (fill_tags->contains(std::string_view(bam_aux_tag(curr_aux), 2))) {
      mAuxTags.emplace_back(AuxTag(curr_aux));
    }

    curr_aux = bam_aux_next(aln, curr_aux);
    if (errno == EINVAL) {
      throw std::runtime_error("aux data for BAM/CRAM record is corrupt");
    }
  }
}

auto Alignment::CigarString() const -> std::string {
  std::string result;
  result.reserve(mCigar.size() * 4);

  for (const auto& unit : mCigar) {
    absl::StrAppend(&result, bam_cigar_oplen(unit), std::string(1, bam_cigar_opchr(unit)));
  }

  result.shrink_to_fit();
  return result;
}

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
         mSamFlag == 0 && mMapQual == 0 && mQname.empty() && mSeq.empty() && mQual.empty() && mCigar.empty() &&
         mAuxTags.empty();
}

auto Alignment::TagNamesView() const -> std::vector<std::string_view> {
  std::vector<std::string_view> result;
  result.reserve(mAuxTags.size());
  std::transform(mAuxTags.cbegin(), mAuxTags.cend(), std::back_inserter(result), std::mem_fn(&AuxTag::Name));
  std::sort(result.begin(), result.end());
  return result;
}

auto Alignment::FindTag(std::string_view tag_name) const -> TagsConstIterator {
  return std::find_if(mAuxTags.cbegin(), mAuxTags.cend(),
                      [&tag_name](const AuxTag& curr_tag) -> bool { return curr_tag.Name() == tag_name; });
}

auto Alignment::ToString(const Reference& ref) const -> std::string {
  const auto chrom = ref.FindChromByIndex(mChromIdx);
  const auto mate_chrom = ref.FindChromByIndex(mMateChromIdx);

  const auto both_chroms_same = chrom.ok() && mate_chrom.ok() && chrom->Name() == mate_chrom->Name();
  const auto rnext = both_chroms_same ? "=" : mate_chrom.ok() ? mate_chrom->Name() : "*";

  std::string fastq_quality;
  fastq_quality.reserve(mQual.size());
  static constexpr u8 PHRED_QUALITY_OFFSET = 33;
  for (const auto& base_qual : mQual) {
    fastq_quality.push_back(static_cast<char>(base_qual + PHRED_QUALITY_OFFSET));
  }

  std::string tags_data;
  static constexpr usize NUM_ESTIMATED_CHARS_PER_TAG = 2048;
  tags_data.reserve(mAuxTags.size() * NUM_ESTIMATED_CHARS_PER_TAG);
  for (const auto& tag : mAuxTags) {
    absl::StrAppend(&tags_data, "\t", tag.ToString());
  }
  tags_data.shrink_to_fit();

  return fmt::format(
      "{QNAME}\t{FLAG}\t{RNAME}\t{POS}\t{MAPQ}\t{CIGAR}\t{RNEXT}\t{PNEXT}\t{TLEN}\t{SEQ}\t{QUAL}{TAGS_IF_PRESENT}\n",
      fmt::arg("QNAME", mQname.empty() ? "*" : mQname), fmt::arg("FLAG", mSamFlag),
      fmt::arg("RNAME", chrom.ok() ? chrom->Name() : "*"), fmt::arg("POS", mStart0 >= 0 ? mStart0 + 1 : 0),
      fmt::arg("MAPQ", mMapQual), fmt::arg("CIGAR", CigarString()), fmt::arg("RNEXT", rnext),
      fmt::arg("PNEXT", mMateStart0 >= 0 ? mMateStart0 + 1 : 0), fmt::arg("TLEN", mInsertSize), fmt::arg("SEQ", mSeq),
      fmt::arg("QUAL", fastq_quality), fmt::arg("TAGS_IF_PRESENT", tags_data));
}

auto Alignment::operator==(const Alignment& rhs) const -> bool {
  return mStart0 == rhs.mStart0 && mMateStart0 == rhs.mMateStart0 && mInsertSize == rhs.mInsertSize &&
         mChromIdx == rhs.mChromIdx && mMateChromIdx == rhs.mMateChromIdx && mSamFlag == rhs.mSamFlag &&
         mMapQual == rhs.mMapQual && mQname == rhs.mQname && mSeq == rhs.mSeq && mQual == rhs.mQual &&
         mCigar == rhs.mCigar && mAuxTags == rhs.mAuxTags;
}

auto Alignment::operator!=(const Alignment& rhs) const -> bool { return !(rhs == *this); }

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto Alignment::GetSoftClips(std::vector<u32>* clip_sizes, std::vector<u32>* read_positions,
                             std::vector<u32>* genome_positions, bool use_padded) const -> bool {
  // initialize positions & flags
  auto ref_position = static_cast<u32>(mStart0);
  u32 read_position = 0;
  bool soft_clip_found = false;
  bool first_cigar_op = true;

  for (const auto& val : mCigar) {
    const auto cig_unit = CigarUnit(val);
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
