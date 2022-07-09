#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "lancet2/cli_params.h"
#include "lancet2/fasta_reader.h"
#include "lancet2/ref_window.h"

namespace lancet2 {
using WindowPtr = std::shared_ptr<RefWindow>;

class WindowBuilder {
 public:
  WindowBuilder(const std::filesystem::path& ref, u32 region_padding, u32 window_length, u32 pct_window_overlap);
  WindowBuilder() = delete;

  [[nodiscard]] auto AddRegion(const std::string& region_str) -> absl::Status;
  [[nodiscard]] auto AddRegionsBatch(absl::Span<const std::string> region_strs) -> absl::Status;
  [[nodiscard]] auto AddRegionsFromBedFile(const std::filesystem::path& bed) -> absl::Status;

  void AddAllRefRegions();

  [[nodiscard]] auto GetSize() const noexcept { return inputRegions.size(); }
  [[nodiscard]] auto IsEmpty() const noexcept { return inputRegions.empty(); }

  /// 1. Combine the input regions from bed file & samtools-style regions.
  ///    Non-OK status is returned if no input regions have been added.
  /// 2. Add `regionPadding` to each input region from the previous step
  /// 3. Build result windows each `windowLength` in length and
  ///    overlap of `pctWindowOverlap`% between consecutive windows
  [[nodiscard]] auto BuildWindows() const -> absl::StatusOr<std::vector<WindowPtr>>;

  [[nodiscard]] static auto StepSize(u32 pct_overlap, u32 window_length) -> i64;

 private:
  FastaReader refRdr;
  u32 regionPadding = DEFAULT_REGION_PAD_LENGTH;
  u32 windowLength = DEFAULT_WINDOW_LENGTH;
  u32 pctWindowOverlap = DEFAULT_PCT_WINDOW_OVERLAP;
  absl::flat_hash_set<std::string> inputRegions;

  // Parse samtools region string, padding is not added yet
  [[nodiscard]] static auto ParseRegion(std::string_view region_str) -> absl::StatusOr<RefWindow>;

  // Parse input regions from bed file, padding is not added yet
  [[nodiscard]] static auto ParseBed(const std::filesystem::path& bed) -> absl::StatusOr<std::vector<RefWindow>>;

  // Add `regionPadding` to start and end, while checking for coordinate under/over-flow
  [[nodiscard]] auto PadWindow(const RefWindow& w) const -> absl::StatusOr<RefWindow>;
};

/// 1. Combine the input regions from bed file & samtools-style regions.
///    Non-OK status is returned if no input regions have been added.
/// 2. Add `regionPadding` to each input region from the previous step
/// 3. Build result windows each `windowLength` in length and
///    overlap of `pctWindowOverlap`% between consecutive windows
/// NOTE: calls std::exit on failure
[[nodiscard]] auto BuildWindows(const absl::flat_hash_map<std::string, i64>& contig_ids, const CliParams& params)
    -> std::vector<WindowPtr>;
}  // namespace lancet2
