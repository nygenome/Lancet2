#ifndef SRC_LANCET_CORE_READ_COLLECTOR_H_
#define SRC_LANCET_CORE_READ_COLLECTOR_H_

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/read.h"
#include "lancet/core/sample_info.h"
#include "lancet/hts/alignment.h"
#include "lancet/hts/extractor.h"
#include "lancet/hts/reference.h"

namespace lancet::core {

class ReadCollector {
 public:
  static constexpr f64 DEFAULT_MAX_WINDOW_COVERAGE = 500.0;

  struct Params {
    // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
    std::filesystem::path mRefPath;
    std::vector<std::filesystem::path> mNormalPaths;
    std::vector<std::filesystem::path> mTumorPaths;

    f64 mMaxSampleCov = DEFAULT_MAX_WINDOW_COVERAGE;
    bool mNoCtgCheck = false;
    bool mExtractPairs = false;
    // NOLINTEND(misc-non-private-member-variables-in-classes)

    [[nodiscard]] auto SamplesCount() const -> usize { return mNormalPaths.size() + mTumorPaths.size(); }
  };

  using Read = cbdg::Read;
  using Region = hts::Reference::Region;

  explicit ReadCollector(Params params);

  struct Result {
    std::vector<Read> mSampleReads;
    std::vector<SampleInfo> mSampleList;
  };

  [[nodiscard]] auto CollectRegionResult(const Region& region) -> Result;
  [[nodiscard]] auto IsGermlineMode() const noexcept -> bool { return mIsGermlineMode; }

  [[nodiscard]] static auto IsActiveRegion(const Params& params, const Region& region) -> bool;
  [[nodiscard]] static auto BuildSampleNameList(const Params& params) -> std::vector<std::string>;

 private:
  using ExtractorPtr = std::unique_ptr<hts::Extractor>;
  using AlnAndRefPaths = std::array<std::filesystem::path, 2>;
  using MateRegionsMap = absl::flat_hash_map<std::string, hts::Alignment::MateInfo>;
  using SampleExtractors = absl::flat_hash_map<SampleInfo, ExtractorPtr, SampleInfo::Hash, SampleInfo::Equal>;

  Params mParams;
  bool mIsGermlineMode;
  SampleExtractors mExtractors;
  std::vector<SampleInfo> mSampleList;

  [[nodiscard]] static auto MakeSampleList(const Params& params) -> std::vector<SampleInfo>;

  using MateNameAndLocation = std::pair<std::string, hts::Alignment::MateInfo>;
  [[nodiscard]] static auto RevSortMateRegions(const MateRegionsMap& data) -> std::vector<MateNameAndLocation>;
  [[nodiscard]] static auto MakeRegSpec(const hts::Alignment::MateInfo& info, const hts::Extractor* ext) -> std::string;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_READ_COLLECTOR_H_
