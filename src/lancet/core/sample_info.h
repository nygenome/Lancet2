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

  [[nodiscard]] auto Path() const noexcept -> std::filesystem::path const& { return mFilePath; }
  [[nodiscard]] auto FileName() const -> std::string { return mFilePath.filename().string(); }

  [[nodiscard]] auto TagKind() const noexcept -> cbdg::Label::Tag { return mTag; }
  [[nodiscard]] auto NumSampledReads() const noexcept -> u64 { return mNumSampledReads; }
  [[nodiscard]] auto NumSampledBases() const noexcept -> u64 { return mNumSampledBases; }
  [[nodiscard]] auto SampleName() const noexcept -> std::string_view { return mSampleName; }

  /// Per-sample mean read depth over a window of length window_length.
  /// Used by SDFC (Site Depth Fold Change) to normalize per-sample site depth.
  [[nodiscard]] auto MeanCoverage(u64 const window_length) const noexcept -> f64 {
    return window_length > 0 ? static_cast<f64>(mNumSampledBases) / static_cast<f64>(window_length)
                             : 0.0;
  }

  /// Cross-sample mean read depth: sum of all samples' bases / window_length.
  [[nodiscard]] static auto CrossSampleMeanCoverage(absl::Span<SampleInfo const> samples,
                                                    u64 const window_length) -> f64 {
    static auto const SUMMER = [](u64 const sum, SampleInfo const& sinfo) -> u64 {
      return sum + sinfo.NumSampledBases();
    };

    u64 const total_bases = std::accumulate(samples.cbegin(), samples.cend(), u64{0}, SUMMER);
    return static_cast<f64>(total_bases) / static_cast<f64>(window_length);
  }

  friend auto operator<(SampleInfo const& lhs, SampleInfo const& rhs) -> bool {
    return lhs.TagKind() != rhs.TagKind()
               ? static_cast<u8>(lhs.TagKind()) < static_cast<u8>(rhs.TagKind())
               : lhs.SampleName() < rhs.SampleName();
  }

  struct Hash {
    using IsTransparent = void;

    auto operator()(SampleInfo const& sinfo) const -> usize {
      return lancet::base::HashStr64(sinfo.SampleName());
    }
    auto operator()(std::string_view sample_name) const -> usize {
      return lancet::base::HashStr64(sample_name);
    }
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

  /// 0-based index assigned after sorting by (TagKind, SampleName).
  /// Multiple BAM/CRAM files sharing the same SM read group tag and same role
  /// receive the same index — they are one logical sample split across files.
  ///
  /// This index determines three things:
  ///   1. Which entry in Node::mCounts tracks this sample's read support
  ///   2. The position of this sample's FORMAT column in VCF output
  ///   3. Which bit is set in the SampleMask for N-sample graph coloring
  ///
  /// Assigned by MakeSampleList() after sorting — never set by the constructor.
  [[nodiscard]] auto SampleIndex() const noexcept -> usize { return mSampleIndex; }

  void SetSampleIndex(usize const index) { mSampleIndex = index; }
  void SetNumSampledReads(u64 const num_reads) { mNumSampledReads = num_reads; }
  void SetNumSampledBases(u64 const num_bases) { mNumSampledBases = num_bases; }

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  u64 mNumSampledReads = 0;
  u64 mNumSampledBases = 0;
  usize mSampleIndex = 0;
  std::string mSampleName;
  std::filesystem::path mFilePath;
  // ── 1B Align ────────────────────────────────────────────────────────────
  cbdg::Label::Tag mTag = cbdg::Label::REFERENCE;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_SAMPLE_INFO_H_
