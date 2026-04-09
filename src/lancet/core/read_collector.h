#ifndef SRC_LANCET_CORE_READ_COLLECTOR_H_
#define SRC_LANCET_CORE_READ_COLLECTOR_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/read.h"
#include "lancet/core/sample_info.h"
#include "lancet/hts/alignment.h"
#include "lancet/hts/extractor.h"
#include "lancet/hts/reference.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace lancet::core {

class ReadCollector {
 public:
  static constexpr f64 DEFAULT_MAX_WINDOW_COVERAGE = 1000.0;

  struct Params {
    // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
    std::filesystem::path mRefPath;
    std::vector<std::filesystem::path> mNormalPaths;
    std::vector<std::filesystem::path> mTumorPaths;

    f64 mMaxSampleCov = DEFAULT_MAX_WINDOW_COVERAGE;
    bool mNoCtgCheck = false;
    bool mExtractPairs = false;
    // NOLINTEND(misc-non-private-member-variables-in-classes)

    [[nodiscard]] auto SamplesCount() const -> usize {
      return mNormalPaths.size() + mTumorPaths.size();
    }
  };

  using Read = cbdg::Read;
  using Region = hts::Reference::Region;

  explicit ReadCollector(Params params);

  struct Result {
    std::vector<Read> mSampleReads;
    std::vector<SampleInfo> mSampleList;
  };

  [[nodiscard]] auto CollectRegionResult(Region const& region) -> Result;
  [[nodiscard]] auto IsTumorNormalMode() const noexcept -> bool { return mIsTumorNormalMode; }

  [[nodiscard]] static auto IsActiveRegion(Params const& params, Region const& region) -> bool;
  [[nodiscard]] static auto BuildSampleNameList(Params const& params) -> std::vector<std::string>;

 private:
  using ExtractorPtr = std::unique_ptr<hts::Extractor>;
  using AlnAndRefPaths = std::array<std::filesystem::path, 2>;
  using SampleExtractors =
      absl::flat_hash_map<SampleInfo, ExtractorPtr, SampleInfo::Hash, SampleInfo::Equal>;

  /// Maps qname hash (u64) -> mate location info for out-of-region mate retrieval.
  /// Using u64 hashes instead of std::string keys avoids string copies during Pass 1.
  using MateRegionsMap = absl::flat_hash_map<u64, hts::Alignment::MateInfo>;

  Params mParams;
  bool mIsTumorNormalMode{false};
  SampleExtractors mExtractors;
  std::vector<SampleInfo> mSampleList;

  [[nodiscard]] static auto MakeSampleList(Params const& params) -> std::vector<SampleInfo>;

  /// Computes a deterministic u64 hash for a query name string_view.
  [[nodiscard]] static auto HashQname(std::string_view qname) -> u64;

  using MateHashAndLocation = std::pair<u64, hts::Alignment::MateInfo>;
  [[nodiscard]] static auto RevSortMateRegions(MateRegionsMap const& data)
      -> std::vector<MateHashAndLocation>;
  [[nodiscard]] static auto MakeRegSpec(hts::Alignment::MateInfo const& info,
                                        hts::Extractor const* ext) -> std::string;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_READ_COLLECTOR_H_
