#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "absl/status/status.h"
#include "lancet/cli_params.h"
#include "lancet/fasta_reader.h"
#include "lancet/ref_window.h"
#include "lancet/statusor.h"

namespace lancet {
using WindowPtr = std::shared_ptr<RefWindow>;

class WindowBuilder {
 public:
  WindowBuilder(const std::filesystem::path& ref, std::uint32_t region_padding, std::uint32_t window_length,
                std::uint32_t pct_window_overlap);
  WindowBuilder() = delete;

  [[nodiscard]] auto AddSamtoolsRegion(const std::string& region_str) -> absl::Status;
  [[nodiscard]] auto AddBedFileRegions(const std::filesystem::path& bed) -> absl::Status;

  void AddAllRefRegions();

  [[nodiscard]] auto Size() const noexcept { return inputRegions.size(); }
  [[nodiscard]] auto IsEmpty() const noexcept { return inputRegions.empty(); }

  /// 1. Combine the input regions from bed file & samtools-style regions.
  ///    Non-OK status is returned if no input regions have been added.
  /// 2. Add `regionPadding` to each input region from the previous step
  /// 3. Build result windows each `windowLength` in length and
  ///    overlap of `pctWindowOverlap`% between consecutive windows
  [[nodiscard]] auto BuildWindows(const absl::flat_hash_map<std::string, std::int64_t>& contig_ids) const
      -> StatusOr<std::vector<WindowPtr>>;

  [[nodiscard]] static auto StepSize(std::uint32_t pct_overlap, std::uint32_t window_length) -> std::int64_t;

 private:
  FastaReader refRdr;
  std::uint32_t regionPadding = DEFAULT_REGION_PAD_LENGTH;
  std::uint32_t windowLength = DEFAULT_WINDOW_LENGTH;
  std::uint32_t pctWindowOverlap = DEFAULT_PCT_WINDOW_OVERLAP;
  std::vector<RefWindow> inputRegions;  // sequence only added in `BuildWindows`

  // Parse samtools region string, padding is not added yet
  [[nodiscard]] static auto ParseRegion(std::string_view region_str) -> StatusOr<RefWindow>;

  // Parse input regions from bed file, padding is not added yet
  [[nodiscard]] static auto ParseBed(const std::filesystem::path& bed) -> StatusOr<std::vector<RefWindow>>;

  // Add `regionPadding` to start and end, while checking for coordinate under/over-flow
  [[nodiscard]] auto PadWindow(const RefWindow& w) const -> StatusOr<RefWindow>;
};

/// 1. Combine the input regions from bed file & samtools-style regions.
///    Non-OK status is returned if no input regions have been added.
/// 2. Add `regionPadding` to each input region from the previous step
/// 3. Build result windows each `windowLength` in length and
///    overlap of `pctWindowOverlap`% between consecutive windows
/// NOTE: calls std::exit on failure
[[nodiscard]] auto BuildWindows(const absl::flat_hash_map<std::string, std::int64_t>& contig_ids,
                                const CliParams& params) -> std::vector<WindowPtr>;
}  // namespace lancet
