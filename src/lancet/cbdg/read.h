#ifndef SRC_LANCET_CBDG_READ_H_
#define SRC_LANCET_CBDG_READ_H_

#include <cmath>
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
        mQuality(aln.QualView().cbegin(), aln.QualView().cend()) {
    static constexpr u8 DEFAULT_MIN_READ_MAP_QUAL = 20;
    if (aln.MapQual() < DEFAULT_MIN_READ_MAP_QUAL) {
      mPassesAlnFilters = false;
    }

    // AS: Alignment score; XS: Suboptimal alignment score
    static constexpr f64 DEFAULT_MIN_READ_AS_XS_PCT_DIFF = 1.0;
    if (aln.HasTag("AS") && aln.HasTag("XS")) {
      const auto as_tag = aln.GetTag<i64>("AS").value();
      const auto xs_tag = aln.GetTag<i64>("XS").value();
      const auto pct_score_decrease = (static_cast<f64>(std::abs(xs_tag - as_tag)) * 100.0) / static_cast<f64>(as_tag);
      mPctAlnScoresDiff = static_cast<u8>(std::round(pct_score_decrease));

      if (pct_score_decrease < DEFAULT_MIN_READ_AS_XS_PCT_DIFF) {
        mPassesAlnFilters = false;
      }
    }

    // XT type: Unique/Repeat/N/Mate-sw
    // XT:A:M (one-mate recovered) means that one of the pairs is uniquely mapped and the other isn't
    // Heng Li: If the read itself is a repeat and can't be mapped without relying on its mate, you
    // see "XT:Z:R". Nonetheless, the mapping quality is not necessarily zero. When its mate can be
    // mapped unambiguously, the read can still be mapped confidently and thus assigned a high mapQ.
    // MapQ is computed for the read pair. XT is determined from a single read.
    // XA -- BWA (Illumina): alternative hits; format: (chr,pos,CIGAR,NM;)
    // return aln.HasTag("XT") || aln.HasTag("XA");
  }

  [[nodiscard]] auto StartPos0() const noexcept -> i64 { return mStart0; }
  [[nodiscard]] auto ChromIndex() const noexcept -> i32 { return mChromIdx; }
  [[nodiscard]] auto BitwiseFlag() const noexcept -> hts::Alignment::BitwiseFlag { return mSamFlag; }
  [[nodiscard]] auto MapQual() const noexcept -> u8 { return mMapQual; }

  [[nodiscard]] auto PassesAlnFilters() const noexcept -> bool { return mPassesAlnFilters; }
  [[nodiscard]] auto PctAlnScoresDiff() const noexcept -> u8 { return mPctAlnScoresDiff; }

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
           lhs.mQname == rhs.mQname && lhs.mSequence == rhs.mSequence && lhs.mQuality == rhs.mQuality &&
           lhs.mPassesAlnFilters == rhs.mPassesAlnFilters;
  }

  friend auto operator!=(const Read& lhs, const Read& rhs) noexcept -> bool { return !(lhs == rhs); }

 private:
  i64 mStart0 = -1;
  i32 mChromIdx = -1;
  u16 mSamFlag = 0;
  u8 mMapQual = 0;

  u8 mPctAlnScoresDiff = 100.0;
  bool mPassesAlnFilters = true;

  Label::Tag mTag;
  std::string mQname;
  std::string mSequence;
  std::string mSampleName;
  std::vector<u8> mQuality;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_READ_H_
