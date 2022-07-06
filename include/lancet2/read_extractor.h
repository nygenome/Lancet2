#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <string_view>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "lancet2/cli_params.h"
#include "lancet2/core_enums.h"
#include "lancet2/genomic_region.h"
#include "lancet2/hts_reader.h"
#include "lancet2/read_info.h"

namespace lancet2 {
using ReadInfoList = std::vector<ReadInfo>;

class ReadExtractor {
 public:
  explicit ReadExtractor(std::shared_ptr<const CliParams> p);
  ReadExtractor() = delete;

  void SetTargetRegion(const GenomicRegion& region);

  [[nodiscard]] auto IsActiveRegion() const -> bool { return isActiveRegion; }
  [[nodiscard]] auto AverageCoverage() const -> double { return avgCoverage; }

  [[nodiscard]] auto Extract() -> ReadInfoList;

 private:
  HtsReader tmrRdr;
  HtsReader nmlRdr;
  GenomicRegion targetRegion{"\0", -1, -1};
  bool isActiveRegion = false;
  double avgCoverage = 0.0F;
  std::shared_ptr<const CliParams> params;

  [[nodiscard]] auto ExtractReads(HtsReader* rdr, ReadInfoList* result, SampleLabel label)
      -> absl::flat_hash_map<std::string, GenomicRegion>;

  void ExtractPairs(HtsReader* rdr, const absl::flat_hash_map<std::string, GenomicRegion>& mate_info,
                    ReadInfoList* result, SampleLabel label);

  [[nodiscard]] static auto PassesFilters(const HtsAlignment& aln, const CliParams& params) -> bool;
  [[nodiscard]] static auto PassesTmrFilters(const HtsAlignment& aln, const CliParams& params) -> bool;

  struct EvalResult {
    double coverage = 0.0F;
    bool isActiveRegion = false;
  };

  [[nodiscard]] static auto EvaluateRegion(HtsReader* rdr, const GenomicRegion& region, const CliParams& params)
      -> EvalResult;

  static void FillMDMismatches(std::string_view md, std::string_view quals, std::int64_t aln_start,
                               std::uint32_t min_bq, std::map<std::uint32_t, std::uint32_t>* result);
};
}  // namespace lancet2
