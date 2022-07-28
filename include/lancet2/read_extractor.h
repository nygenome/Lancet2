#pragma once

#include <map>
#include <memory>
#include <string_view>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "lancet2/cli_params.h"
#include "lancet2/core_enums.h"
#include "lancet2/fractional_sampler.h"
#include "lancet2/genomic_region.h"
#include "lancet2/hts_reader.h"
#include "lancet2/read_info.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
using ReadInfoList = std::vector<ReadInfo>;

class ReadExtractor {
 public:
  explicit ReadExtractor(std::shared_ptr<const CliParams> p);
  ReadExtractor() = delete;

  struct ScanRegionResult {
    u64 NumReadBases = 0;
    double AverageCoverage = 0.0;
    bool HasMutationEvidence = false;
  };

  [[nodiscard]] auto ScanRegion(const GenomicRegion& region) -> ScanRegionResult;
  [[nodiscard]] auto ExtractReads(const GenomicRegion& region, double sampleFraction = 1.0) -> ReadInfoList;

 private:
  std::vector<std::unique_ptr<HtsReader>> readers;
  std::vector<SampleLabel> labels;
  std::shared_ptr<const CliParams> params;
  FractionalSampler sampler;

  using MateInfoMap = absl::flat_hash_map<std::string, GenomicRegion>;
  [[nodiscard]] auto FetchReads(usize sampleIdx, const GenomicRegion& region, ReadInfoList* result) -> MateInfoMap;
  void FetchPairs(usize sampleIdx, const GenomicRegion& region, const MateInfoMap& mate_info, ReadInfoList* result);

  [[nodiscard]] auto ScanSampleRegion(usize sampleIdx, const GenomicRegion& region) -> ScanRegionResult;

  [[nodiscard]] static auto PassesFilters(const HtsAlignment& aln, const CliParams& params, SampleLabel label) -> bool;

  static void FillMD(std::string_view md, std::string_view quals, i64 aln_start, u32 min_bq,
                     std::map<u32, u32>* result);
};
}  // namespace lancet2
