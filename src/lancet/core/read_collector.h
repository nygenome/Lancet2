#ifndef SRC_LANCET_CORE_READ_COLLECTOR_H_
#define SRC_LANCET_CORE_READ_COLLECTOR_H_

#include <filesystem>
#include <memory>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"
#include "lancet/base/downsampler.h"
#include "lancet/base/hash.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/read.h"
#include "lancet/core/sample_info.h"
#include "lancet/hts/extractor.h"
#include "lancet/hts/reference.h"

namespace lancet::core {

class ReadCollector {
 public:
  static constexpr f64 DEFAULT_MAX_WINDOW_COVERAGE = 1000.0;

  using InsertRange = std::array<i64, 2>;
  using BamCramWithInsert = std::pair<std::filesystem::path, InsertRange>;

  struct Params {
    // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
    std::filesystem::path mRefPath;
    std::vector<BamCramWithInsert> mNormals;
    std::vector<BamCramWithInsert> mTumors;

    f64 mMaxWinCov = DEFAULT_MAX_WINDOW_COVERAGE;
    bool mNoCtgCheck = false;
    bool mExtractPairs = false;
    // NOLINTEND(misc-non-private-member-variables-in-classes)

    [[nodiscard]] auto SamplesCount() const -> usize { return mNormals.size() + mTumors.size(); }
  };

  using Read = cbdg::Read;
  using Region = hts::Reference::Region;

  explicit ReadCollector(Params params);

  struct Result {
    std::vector<Read> mSampleReads;
    std::vector<SampleInfo> mSampleList;
  };

  [[nodiscard]] auto CollectRegionResult(const Region& region) -> Result;

  using AlnAndRefPaths = std::array<std::filesystem::path, 2>;
  [[nodiscard]] static auto EstimateInsertRange(const AlnAndRefPaths& paths) -> InsertRange;
  [[nodiscard]] static auto IsActiveRegion(const Params& params, const Region& region) -> bool;
  [[nodiscard]] static auto BuildSampleNameList(const Params& params) -> std::vector<std::string>;

 private:
  using ExtractorPtr = std::unique_ptr<hts::Extractor>;
  using MateRegionsMap = absl::flat_hash_map<std::string, hts::Alignment::MateInfo>;
  using SampleExtractors = absl::flat_hash_map<SampleInfo, ExtractorPtr, SampleInfo::Hash, SampleInfo::Equal>;

  Params mParams;
  bool mIsGermlineMode;
  Downsampler mDownsampler;
  SampleExtractors mExtractors;
  std::vector<SampleInfo> mSampleList;

  [[nodiscard]] auto EstimateCoverage(const SampleInfo& sinfo, const Region& region) const -> f64;

  [[nodiscard]] static auto FailsTier1Check(const hts::Alignment& aln) -> bool;
  [[nodiscard]] static auto FailsTier2Check(const hts::Alignment& aln) -> bool;
  [[nodiscard]] static auto MakeSampleList(const Params& params) -> std::vector<SampleInfo>;

  [[nodiscard]] static auto BuildSortedMateInfos(const MateRegionsMap& data) -> std::vector<hts::Alignment::MateInfo>;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_READ_COLLECTOR_H_
