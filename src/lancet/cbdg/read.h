#ifndef SRC_LANCET_CBDG_READ_H_
#define SRC_LANCET_CBDG_READ_H_

#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/label.h"
#include "lancet/hts/alignment.h"

namespace lancet::cbdg {

class Read {
 public:
  explicit Read(const hts::Alignment& aln, std::string sample_name, const Label::Tag tag)
      : mStart0(aln.StartPos0()), mChromIdx(aln.ChromIndex()), mSamFlag(aln.FlagRaw()), mMapQual(aln.MapQual()),
        mTag(tag), mQname(aln.QnameView()), mSequence(aln.SeqView()), mSampleName(std::move(sample_name)),
        mQuality(aln.QualView().cbegin(), aln.QualView().cend()) {}

  [[nodiscard]] auto StartPos0() const noexcept -> i64 { return mStart0; }
  [[nodiscard]] auto ChromIndex() const noexcept -> i32 { return mChromIdx; }
  [[nodiscard]] auto BitwiseFlag() const noexcept -> hts::Alignment::BitwiseFlag { return mSamFlag; }
  [[nodiscard]] auto MapQual() const noexcept -> u8 { return mMapQual; }

  [[nodiscard]] auto SrcLabel() const noexcept -> Label { return {mTag}; }
  [[nodiscard]] auto TagKind() const noexcept -> Label::Tag { return mTag; }
  [[nodiscard]] auto QnamePtr() const noexcept -> const char* { return mQname.c_str(); }
  [[nodiscard]] auto SeqPtr() const noexcept -> const char* { return mSequence.c_str(); }
  [[nodiscard]] auto QnameView() const noexcept -> std::string_view { return mQname; }
  [[nodiscard]] auto SeqView() const noexcept -> std::string_view { return mSequence; }
  [[nodiscard]] auto QualView() const noexcept -> absl::Span<const u8> { return mQuality; }
  [[nodiscard]] auto Length() const noexcept -> usize { return mSequence.size(); }
  [[nodiscard]] auto SampleName() const noexcept -> std::string_view { return mSampleName; }

  template <typename HashState>
  friend auto AbslHashValue(HashState hash_state, const Read& read) -> HashState {
    return HashState::combine(std::move(hash_state), read.mSampleName, read.mStart0, read.mChromIdx, read.mSamFlag,
                              read.mMapQual, static_cast<u8>(read.mTag), read.mQname, read.mSequence, read.mQuality);
  }

  friend auto operator==(const Read& lhs, const Read& rhs) noexcept -> bool {
    return lhs.mSampleName == rhs.mSampleName && lhs.mStart0 == rhs.mStart0 && lhs.mChromIdx == rhs.mChromIdx &&
           lhs.mSamFlag == rhs.mSamFlag && lhs.mMapQual == rhs.mMapQual && lhs.mTag == rhs.mTag &&
           lhs.mQname == rhs.mQname && lhs.mSequence == rhs.mSequence && lhs.mQuality == rhs.mQuality;
  }

  friend auto operator!=(const Read& lhs, const Read& rhs) noexcept -> bool { return !(lhs == rhs); }

 private:
  i64 mStart0 = -1;
  i32 mChromIdx = -1;
  u16 mSamFlag = 0;
  u8 mMapQual = 0;

  Label::Tag mTag;
  std::string mQname;
  std::string mSequence;
  std::string mSampleName;
  std::vector<u8> mQuality;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_READ_H_
