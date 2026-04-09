#ifndef SRC_LANCET_CORE_SAMPLE_INFO_H_
#define SRC_LANCET_CORE_SAMPLE_INFO_H_

#include "lancet/base/hash.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/label.h"

#include "absl/types/span.h"

#include <filesystem>
#include <numeric>
#include <string>
#include <string_view>
#include <utility>

namespace lancet::core {

class SampleInfo {
 public:
  explicit SampleInfo(std::string sample_name, std::filesystem::path fpath,
                      cbdg::Label::Tag const tag)
      : mSampleName(std::move(sample_name)), mFilePath(std::move(fpath)), mTag(tag) {}

  [[nodiscard]] auto Path() const noexcept -> std::filesystem::path { return mFilePath; }
  [[nodiscard]] auto FileName() const noexcept -> std::string {
    return mFilePath.filename().string();
  }
  [[nodiscard]] auto TagKind() const noexcept -> cbdg::Label::Tag { return mTag; }

  [[nodiscard]] auto NumSampledReads() const noexcept -> u64 { return mNumSampledReads; }
  [[nodiscard]] auto NumSampledBases() const noexcept -> u64 { return mNumSampledBases; }
  [[nodiscard]] auto SampleName() const noexcept -> std::string_view { return mSampleName; }

  [[nodiscard]] static auto CombinedSampledCov(absl::Span<SampleInfo const> samples,
                                               u64 const ref_len) -> f64 {
    static auto const SUMMER = [](u64 const sum, SampleInfo const& sinfo) -> u64 {
      return sum + sinfo.NumSampledBases();
    };

    u64 const total_bases = std::accumulate(samples.cbegin(), samples.cend(), 0, SUMMER);
    return static_cast<f64>(total_bases) / static_cast<f64>(ref_len);
  }

  friend auto operator<(SampleInfo const& lhs, SampleInfo const& rhs) -> bool {
    return lhs.TagKind() != rhs.TagKind()
               ? static_cast<u8>(lhs.TagKind()) < static_cast<u8>(rhs.TagKind())
               : lhs.SampleName() < rhs.SampleName();
  }

  struct Hash {
    using IsTransparent = void;

    auto operator()(SampleInfo const& sinfo) const -> usize {
      return HashStr64(sinfo.SampleName());
    }
    auto operator()(std::string_view sample_name) const -> usize { return HashStr64(sample_name); }
  };

  struct Equal {
    using IsTransparent = void;

    auto operator()(SampleInfo const& lhs, SampleInfo const& rhs) const -> bool {
      return lhs.SampleName() == rhs.SampleName();
    }

    auto operator()(SampleInfo const& lhs, std::string_view rhs_sample) const -> bool {
      return lhs.SampleName() == rhs_sample;
    }

    auto operator()(std::string_view lhs_sample, SampleInfo const& rhs) const -> bool {
      return lhs_sample == rhs.SampleName();
    }
  };

 private:
  u64 mNumSampledReads = 0;
  u64 mNumSampledBases = 0;

  std::string mSampleName;
  std::filesystem::path mFilePath;
  cbdg::Label::Tag mTag = cbdg::Label::REFERENCE;

  friend class ReadCollector;
  void SetNumSampledReads(u64 const num_reads) { mNumSampledReads = num_reads; }
  void SetNumSampledBases(u64 const num_bases) { mNumSampledBases = num_bases; }
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_SAMPLE_INFO_H_
