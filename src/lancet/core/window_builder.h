#ifndef SRC_LANCET_CORE_WINDOW_BUILDER_H_
#define SRC_LANCET_CORE_WINDOW_BUILDER_H_

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/core/window.h"
#include "lancet/hts/reference.h"

namespace lancet::core {

class WindowBuilder {
 public:
  static constexpr u32 DEFAULT_PCT_OVERLAP = 20;
  static constexpr u32 DEFAULT_WINDOW_LENGTH = 1000;
  static constexpr u32 DEFAULT_REGION_PADDING = 500;

  static constexpr u32 MIN_ALLOWED_PCT_OVERLAP = 10;
  static constexpr u32 MAX_ALLOWED_PCT_OVERLAP = 90;
  static constexpr u32 MIN_ALLOWED_WINDOW_LEN = 250;
  static constexpr u32 MAX_ALLOWED_WINDOW_LEN = 2500;
  static constexpr u32 MAX_ALLOWED_REGION_PAD = 1000;

  /// Batch size for pipelined window generation. Used by both WindowBuilder
  /// and PipelineRunner to control memory footprint during WGS runs.
  /// Small enough to avoid OOM, large enough to keep worker threads saturated.
  static constexpr usize BATCH_SIZE = 524288;

  struct Params {
    u32 mWindowLength = DEFAULT_WINDOW_LENGTH;
    u32 mRegionPadding = DEFAULT_REGION_PADDING;
    u32 mPercentOverlap = DEFAULT_PCT_OVERLAP;
  };

  WindowBuilder() = delete;
  explicit WindowBuilder(const std::filesystem::path& ref_path, const Params& params);

  void AddAllReferenceRegions();
  void AddRegion(const std::string& region_spec);
  void AddBatchRegions(absl::Span<const std::string> region_specs);
  void AddBatchRegions(const std::filesystem::path& bed_file);

  [[nodiscard]] auto Size() const noexcept -> usize { return mInputRegions.size(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return mInputRegions.empty(); }

  [[nodiscard]] static auto StepSize(const Params& params) -> i64;

  /// Returns the mathematically computed total number of windows that will be generated
  /// from the current input regions, without performing any allocations.
  [[nodiscard]] auto ExpectedTargetWindows() const -> usize;

  /// Returns all windows at once. Suitable for small region sets (targeted panels).
  [[nodiscard]] auto BuildWindows() const -> std::vector<WindowPtr>;

  /// Returns the next batch of up to BATCH_SIZE windows, starting from the given offset.
  /// The offset parameter tracks progress and should be initialized to 0 on the first call.
  /// Returns an empty vector when all windows have been emitted.
  /// Input regions must be sorted (via SortInputRegions) before the first call.
  [[nodiscard]] auto BuildWindowsBatch(usize& offset) const -> std::vector<WindowPtr>;

  /// Sorts input regions by chromosome index and start position. Must be called before
  /// using BuildWindowsBatch() to ensure sequential emission order.
  void SortInputRegions();

 private:
  Params mParams;
  std::unique_ptr<hts::Reference> mRefPtr;

  using ParseRegionResult = hts::Reference::ParseRegionResult;
  /// Stored as a sorted vector (after SortInputRegions) instead of flat_hash_set,
  /// enabling deterministic sequential batch emission without global post-sort.
  std::vector<ParseRegionResult> mInputRegions;

  // Add `regionPadding` to start and end, while checking for coordinate under/over-flow
  void PadInputRegion(ParseRegionResult& result) const;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_WINDOW_BUILDER_H_
