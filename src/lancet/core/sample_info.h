#ifndef SRC_LANCET_CORE_SAMPLE_INFO_H_
#define SRC_LANCET_CORE_SAMPLE_INFO_H_

#include <algorithm>
#include <filesystem>
#include <numeric>
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

  [[nodiscard]] auto NumSampledReads() const noexcept -> u64 { return mNumSampledReads; }
  [[nodiscard]] auto NumSampledBases() const noexcept -> u64 { return mNumSampledBases; }
  [[nodiscard]] auto MeanTotalCov() const noexcept -> f64 { return mMeanTotalCov; }
  [[nodiscard]] auto MeanSampledCov() const noexcept -> f64 { return mMeanSampledCov; }
  [[nodiscard]] auto PassReadsFraction() const noexcept -> f64 { return mPassReadsFraction; }
  [[nodiscard]] auto SampleName() const noexcept -> std::string_view { return mSampleName; }

  [[nodiscard]] static auto CombinedSampledCov(absl::Span<const SampleInfo> samples, const u64 ref_len) -> f64 {
    static const auto summer = [](const u64 sum, const SampleInfo& sinfo) -> u64 {
      return sum + sinfo.NumSampledBases();
    };
    
    const u64 total_bases = std::accumulate(samples.cbegin(), samples.cend(), 0, summer);
    return static_cast<f64>(total_bases) / static_cast<f64>(ref_len);
  }

  friend auto operator<(const SampleInfo& lhs, const SampleInfo& rhs) -> bool {
    return lhs.TagKind() != rhs.TagKind() ? static_cast<u8>(lhs.TagKind()) < static_cast<u8>(rhs.TagKind())
                                          : lhs.SampleName() < rhs.SampleName();
  }

  struct Hash {
    using is_transparent = void;

    auto operator()(const SampleInfo& sinfo) const -> usize { return HashStr64(sinfo.SampleName()); }
    auto operator()(std::string_view sample_name) const -> usize { return HashStr64(sample_name); }
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
  u64 mNumSampledReads = 0;
  u64 mNumSampledBases = 0;
  f64 mMeanTotalCov = 0.0;
  f64 mMeanSampledCov = 0.0;
  f64 mPassReadsFraction = 0.0;

  std::string mSampleName;
  std::filesystem::path mFilePath;
  cbdg::Label::Tag mTag = cbdg::Label::REFERENCE;

  friend class ReadCollector;
  void SetNumSampledReads(const u64 num_reads) { mNumSampledReads = num_reads; }
  void SetNumSampledBases(const u64 num_bases) { mNumSampledBases = num_bases; }

  void CalculateMeanTotalCov(const u64 total_bases, const u64 ref_len) {
    mMeanTotalCov = static_cast<f64>(total_bases) / static_cast<f64>(ref_len);
  }

  void CalculateMeanSampledCov(const u64 ref_len) {
    mMeanSampledCov = static_cast<f64>(mNumSampledBases) / static_cast<f64>(ref_len);
  }

  void CalculatePassReadsFraction(const u64 pass, const u64 total) {
    mPassReadsFraction = static_cast<f64>(pass) / static_cast<f64>(total);
  }
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_SAMPLE_INFO_H_
