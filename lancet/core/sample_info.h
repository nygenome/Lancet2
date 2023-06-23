#ifndef SRC_LANCET_CORE_SAMPLE_INFO_H_
#define SRC_LANCET_CORE_SAMPLE_INFO_H_

#include <algorithm>
#include <filesystem>
#include <string>
#include <string_view>
#include <utility>

#include "absl/types/span.h"
#include "lancet/base/hash.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/label.h"

namespace lancet::core {

class SampleInfo {
 public:
  explicit SampleInfo(std::string sample_name, std::filesystem::path fpath, const cbdg::Label::Tag tag)
      : mSampleName(std::move(sample_name)), mFilePath(std::move(fpath)), mTag(tag) {}

  [[nodiscard]] auto Path() const noexcept -> std::filesystem::path { return mFilePath; }
  [[nodiscard]] auto FileName() const noexcept -> std::string { return mFilePath.filename().string(); }
  [[nodiscard]] auto TagKind() const noexcept -> cbdg::Label::Tag { return mTag; }

  [[nodiscard]] auto NumReads() const noexcept -> u64 { return mNumReads; }
  [[nodiscard]] auto NumBases() const noexcept -> u64 { return mNumBases; }
  [[nodiscard]] auto MeanCov() const noexcept -> f64 { return mMeanCov; }
  [[nodiscard]] auto SampleName() const noexcept -> std::string_view { return mSampleName; }

  [[nodiscard]] static auto TotalMeanCov(absl::Span<const SampleInfo> samples, const u64 ref_len) -> f64 {
    static const auto summer = [](const u64 sum, const SampleInfo& sinfo) -> u64 { return sum + sinfo.NumBases(); };
    const u64 total_bases = std::accumulate(samples.cbegin(), samples.cend(), 0, summer);
    return static_cast<f64>(total_bases) / static_cast<f64>(ref_len);
  }

  friend auto operator<(const SampleInfo& lhs, const SampleInfo& rhs) -> bool {
    return lhs.TagKind() != rhs.TagKind() ? static_cast<u8>(lhs.TagKind()) < static_cast<u8>(rhs.TagKind())
                                          : lhs.SampleName() < rhs.SampleName();
  }

  struct Hash {
    using is_transparent = void;

    auto operator()(const SampleInfo& sinfo) const -> usize { return HashStr(sinfo.SampleName()); }
    auto operator()(std::string_view sample_name) const -> usize { return HashStr(sample_name); }
  };

  struct Equal {
    using is_transparent = void;

    auto operator()(const SampleInfo& lhs, const SampleInfo& rhs) const -> bool {
      return lhs.SampleName() == rhs.SampleName();
    }

    auto operator()(const SampleInfo& lhs, std::string_view rhs_sample) const -> bool {
      return lhs.SampleName() == rhs_sample;
    }

    auto operator()(std::string_view lhs_sample, const SampleInfo& rhs) const -> bool {
      return lhs_sample == rhs.SampleName();
    }
  };

 private:
  u64 mNumReads = 0;
  u64 mNumBases = 0;
  f64 mMeanCov = 0.0;
  
  i64 mMinExpectedInsert = 0;
  i64 mMaxExpectedInsert = 0;

  std::string mSampleName;
  std::filesystem::path mFilePath;
  cbdg::Label::Tag mTag = cbdg::Label::REFERENCE;

  friend class ReadCollector;
  void SetNumReads(const u64 num_reads) { mNumReads = num_reads; }
  void SetNumBases(const u64 num_bases) { mNumBases = num_bases; }
  void CalculateMeanCov(const u64 ref_len) { mMeanCov = static_cast<f64>(mNumBases) / static_cast<f64>(ref_len); }
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_SAMPLE_INFO_H_
