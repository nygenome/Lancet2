#ifndef SRC_LANCET_HTS_ALIGNMENT_H_
#define SRC_LANCET_HTS_ALIGNMENT_H_

#include <string>
#include <string_view>
#include <utility>
#include <vector>

extern "C" {
#include "htslib/hts.h"
#include "htslib/sam.h"
}

#include "absl/container/flat_hash_set.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/hts/aux_tag.h"
#include "lancet/hts/cigar_unit.h"
#include "lancet/hts/reference.h"

namespace lancet::hts {

class Alignment {
 public:
  enum class Fields : u16 {
    CORE_QNAME = SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_RNEXT | SAM_PNEXT | SAM_TLEN,
    SEQ_QUAL = CORE_QNAME | SAM_SEQ | SAM_QUAL,
    CIGAR_SEQ_QUAL = SEQ_QUAL | SAM_CIGAR,
    AUX_RGAUX = CIGAR_SEQ_QUAL | SAM_AUX | SAM_RGAUX,
  };

  enum class Strand : bool { FWD = true, REV = false };

  class BitwiseFlag {
   public:
    BitwiseFlag(u16 flags) : mFlag(flags) {}
    BitwiseFlag() = default;

    [[nodiscard]] auto GetStrand() const noexcept -> Strand;
    [[nodiscard]] auto GetMateStrand() const noexcept -> Strand;
    [[nodiscard]] auto IsFwdStrand() const noexcept -> bool;
    [[nodiscard]] auto IsRevStrand() const noexcept -> bool;
    [[nodiscard]] auto IsMateFwdStrand() const noexcept -> bool;
    [[nodiscard]] auto IsMateRevStrand() const noexcept -> bool;
    [[nodiscard]] auto IsQcFail() const noexcept -> bool;
    [[nodiscard]] auto IsDuplicate() const noexcept -> bool;
    [[nodiscard]] auto IsPrimary() const noexcept -> bool;
    [[nodiscard]] auto IsSecondary() const noexcept -> bool;
    [[nodiscard]] auto IsSupplementary() const noexcept -> bool;
    [[nodiscard]] auto IsMapped() const noexcept -> bool;
    [[nodiscard]] auto IsUnmapped() const noexcept -> bool;
    [[nodiscard]] auto IsMateMapped() const noexcept -> bool;
    [[nodiscard]] auto IsMateUnmapped() const noexcept -> bool;
    [[nodiscard]] auto IsPairedInSequencing() const noexcept -> bool;
    [[nodiscard]] auto IsMappedProperPair() const noexcept -> bool;
    [[nodiscard]] auto IsRead1() const noexcept -> bool;
    [[nodiscard]] auto IsRead2() const noexcept -> bool;
    [[nodiscard]] auto HasFlagsSet(u16 check_flags) const noexcept -> bool;
    [[nodiscard]] auto HasFlagsUnset(u16 check_flags) const noexcept -> bool;

    [[nodiscard]] operator u16() const noexcept { return mFlag; }

    auto operator==(const BitwiseFlag& rhs) const -> bool { return mFlag == rhs.mFlag; }
    auto operator!=(const BitwiseFlag& rhs) const -> bool { return !(rhs == *this); }

   private:
    u16 mFlag = 0;
  };

  [[nodiscard]] auto StartPos0() const noexcept -> i64 { return mStart0; }
  [[nodiscard]] auto MateStartPos0() const noexcept -> i64 { return mMateStart0; }
  [[nodiscard]] auto InsertSize() const noexcept -> i64 { return mInsertSize; }
  [[nodiscard]] auto ChromIndex() const noexcept -> i32 { return mChromIdx; }
  [[nodiscard]] auto MateChromIndex() const noexcept -> i32 { return mMateChromIdx; }
  [[nodiscard]] auto Flag() const noexcept -> BitwiseFlag { return mSamFlag; }
  [[nodiscard]] auto FlagRaw() const noexcept -> u16 { return mSamFlag; }
  [[nodiscard]] auto MapQual() const noexcept -> u8 { return mMapQual; }

  [[nodiscard]] auto QnameView() const noexcept -> std::string_view { return mQname; }
  [[nodiscard]] auto SeqView() const noexcept -> std::string_view { return mSeq; }
  [[nodiscard]] auto QualView() const noexcept -> absl::Span<const u8> { return mQual; }

  [[nodiscard]] auto CigarData() const -> std::vector<CigarUnit> { return {mCigar.cbegin(), mCigar.cend()}; }
  [[nodiscard]] auto CigarString() const -> std::string;

  struct MateInfo {
    i32 mChromIndex = -1;
    i64 mMateStartPos0 = -1;
  };

  [[nodiscard]] auto MateLocation() const noexcept -> MateInfo;
  [[nodiscard]] auto MateOverlapsRegion(const Reference::Region& region) const noexcept -> bool;

  [[nodiscard]] auto OverlapsRegion(const Reference::Region& region) const noexcept -> bool;

  [[nodiscard]] auto Length() const noexcept -> usize { return mSeq.size(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool;

  [[nodiscard]] auto NumTags() const -> usize { return mAuxTags.size(); }
  [[nodiscard]] auto TagNamesView() const -> std::vector<std::string_view>;

  using TagsConstIterator = std::vector<AuxTag>::const_iterator;
  [[nodiscard]] auto FindTag(std::string_view tag_name) const -> TagsConstIterator;
  [[nodiscard]] auto HasTag(std::string_view tag_name) const -> bool { return FindTag(tag_name) != mAuxTags.cend(); }

  template <typename TagResultValue>
  [[nodiscard]] auto GetTag(std::string_view tag_name) const -> absl::StatusOr<TagResultValue> {
    const auto& itr = FindTag(tag_name);
    if (itr != mAuxTags.cend()) {
      return itr->Value<TagResultValue>();
    }

    const auto msg = fmt::format("Tag {} is not present in the alignment record", tag_name);
    return absl::Status(absl::StatusCode::kNotFound, msg);
  }

  [[nodiscard]] auto GetSoftClips(std::vector<u32>* clip_sizes, std::vector<u32>* read_positions,
                                  std::vector<u32>* genome_positions, bool use_padded = false) const -> bool;

  [[nodiscard]] auto ToString(const Reference& ref) const -> std::string;

  template <typename HashState>
  friend auto AbslHashValue(HashState state, const Alignment& aln) -> HashState {
    return HashState::combine(std::move(state), aln.mStart0, aln.mMateStart0, aln.mInsertSize, aln.mChromIdx,
                              aln.mMateChromIdx, aln.mSamFlag, aln.mMapQual, aln.mQname, aln.mSeq, aln.mQual,
                              aln.mCigar, aln.mAuxTags);
  }

  auto operator==(const Alignment& rhs) const -> bool;
  auto operator!=(const Alignment& rhs) const -> bool;

 private:
  i64 mStart0 = -1;
  i64 mMateStart0 = -1;
  i64 mInsertSize = -1;
  i32 mChromIdx = -1;
  i32 mMateChromIdx = -1;
  u16 mSamFlag = 0;
  u8 mMapQual = 0;
  std::string mQname;

  std::string mSeq;
  std::vector<u8> mQual;
  std::vector<u32> mCigar;
  std::vector<AuxTag> mAuxTags;

  friend class Iterator;
  using TagNamesSet = absl::flat_hash_set<std::string>;

  Alignment() = default;

  void ClearAllFields();
  void PopulateRequestedFields(bam1_t* aln, Alignment::Fields fields, const TagNamesSet* fill_tags);

  void PopulateCoreQname(bam1_t* aln);
  void PopulateCigarSeqQual(bam1_t* aln);
  void PopulateAuxRgAux(bam1_t* aln, const TagNamesSet* fill_tags);
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_ALIGNMENT_H_
