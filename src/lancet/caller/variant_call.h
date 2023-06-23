#ifndef SRC_LANCET_CALLER_VARIANT_CALL_H_
#define SRC_LANCET_CALLER_VARIANT_CALL_H_

#include <array>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_support.h"
#include "lancet/core/sample_info.h"

namespace lancet::caller {

using VariantID = u64;

class VariantCall {
 public:
  static constexpr f64 DEFAULT_MIN_TUMOR_VAF = 0.02;
  static constexpr f64 DEFAULT_MAX_NORMAL_VAF = 0.02;
  static constexpr f64 DEFAULT_MIN_TUMOR_COV = 10.0;
  static constexpr f64 DEFAULT_MIN_NORMAL_COV = 10.0;

  static constexpr u32 DEFAULT_MIN_PHRED_SCORE = 20;
  static constexpr u32 DEFAULT_MIN_TUMOR_ALT_COUNT = 2;
  static constexpr u32 DEFAULT_MAX_NORMAL_ALT_COUNT = 2;

  struct Params {
    f64 mMinTmrVaf = DEFAULT_MIN_TUMOR_VAF;
    f64 mMaxNmlVaf = DEFAULT_MAX_NORMAL_VAF;
    f64 mMinTmrCov = DEFAULT_MIN_TUMOR_COV;
    f64 mMinNmlCov = DEFAULT_MIN_NORMAL_COV;

    u32 mMinPhredScore = DEFAULT_MIN_PHRED_SCORE;
    u32 mMinTmrAltCnt = DEFAULT_MIN_TUMOR_ALT_COUNT;
    u32 mMaxNmlAltCnt = DEFAULT_MAX_NORMAL_ALT_COUNT;
  };

  using Samples = absl::Span<const core::SampleInfo>;
  using Supports = absl::flat_hash_map<std::string_view, std::unique_ptr<VariantSupport>>;
  VariantCall(const RawVariant* var, Supports&& supprts, Samples samps, const Params& prms, usize klen);

  [[nodiscard]] auto ChromIndex() const -> usize { return mChromIndex; }
  [[nodiscard]] auto ChromName() const -> std::string_view { return mChromName; }
  [[nodiscard]] auto StartPos1() const -> usize { return mStartPos1; }
  [[nodiscard]] auto RefAllele() const -> std::string_view { return mRefAllele; }
  [[nodiscard]] auto AltAllele() const -> std::string_view { return mAltAllele; }
  [[nodiscard]] auto Length() const -> i64 { return mVarLength; }
  [[nodiscard]] auto Quality() const -> u8 { return mSiteQuality; }
  [[nodiscard]] auto State() const -> RawVariant::State { return mState; }
  [[nodiscard]] auto Category() const -> RawVariant::Type { return mCategory; }

  [[nodiscard]] auto NumSamples() const -> usize { return mFormatFields.empty() ? 0 : mFormatFields.size() - 1; }
  [[nodiscard]] auto Identifier() const -> VariantID { return mVariantId; }
  [[nodiscard]] auto TotalCoverage() const -> usize { return mTotalSampleCov; }

  [[nodiscard]] auto AsVcfRecord() const -> std::string;

  friend auto operator==(const VariantCall& lhs, const VariantCall& rhs) -> bool {
    return lhs.mVariantId == rhs.mVariantId;
  }
  friend auto operator<(const VariantCall& lhs, const VariantCall& rhs) -> bool {
    // NOLINTBEGIN(readability-braces-around-statements)
    if (lhs.mChromIndex != rhs.mChromIndex) return lhs.mChromIndex < rhs.mChromIndex;
    if (lhs.mStartPos1 != rhs.mStartPos1) return lhs.mStartPos1 < rhs.mStartPos1;
    if (lhs.mRefAllele != rhs.mRefAllele) return lhs.mRefAllele < rhs.mRefAllele;
    if (lhs.mAltAllele != rhs.mAltAllele) return lhs.mAltAllele < rhs.mAltAllele;
    if (lhs.mVarLength != rhs.mVarLength) return lhs.mVarLength < rhs.mVarLength;
    return static_cast<i8>(lhs.mCategory) < static_cast<i8>(rhs.mCategory);
    // NOLINTEND(readability-braces-around-statements)
  }

 private:
  u64 mVariantId;
  usize mChromIndex;
  usize mStartPos1;
  usize mTotalSampleCov;
  std::string mChromName;
  std::string mRefAllele;
  std::string mAltAllele;

  i64 mVarLength;
  u8 mSiteQuality;
  RawVariant::State mState;
  RawVariant::Type mCategory;

  std::string mFilterField;
  std::string mInfoField;
  std::vector<std::string> mFormatFields;

  static constexpr std::string_view REF_HOM = "0/0";
  static constexpr std::string_view HET_ALT = "0/1";
  static constexpr std::string_view ALT_HOM = "1/1";
  static constexpr auto POSSIBLE_GENOTYPES = std::array<std::string_view, 3>{REF_HOM, HET_ALT, ALT_HOM};

  using PerSampleEvidence = absl::flat_hash_map<const core::SampleInfo, std::unique_ptr<VariantSupport>,
                                                core::SampleInfo::Hash, core::SampleInfo::Equal>;

  [[nodiscard]] static auto SomaticScore(const core::SampleInfo& tumor, const PerSampleEvidence& supports) -> u8;
  [[nodiscard]] static auto FirstAndSecondSmallestIndices(const std::array<int, 3>& pls) -> std::array<usize, 2>;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_CALL_H_