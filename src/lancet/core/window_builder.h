#ifndef SRC_LANCET_CORE_WINDOW_BUILDER_H_
#define SRC_LANCET_CORE_WINDOW_BUILDER_H_

#include <filesystem>
#include <memory>
#include <string>

#include "absl/container/flat_hash_set.h"
#include "lancet/base/types.h"
#include "lancet/core/window.h"
#include "lancet/hts/reference.h"

namespace lancet::core {

class WindowBuilder {
 public:
  static constexpr u32 DEFAULT_PCT_OVERLAP = 50;
  static constexpr u32 DEFAULT_WINDOW_LENGTH = 1000;
  static constexpr u32 DEFAULT_REGION_PADDING = 500;

  static constexpr u32 MIN_ALLOWED_PCT_OVERLAP = 50;
  static constexpr u32 MAX_ALLOWED_PCT_OVERLAP = 90;
  static constexpr u32 MIN_ALLOWED_WINDOW_LEN = 500;
  static constexpr u32 MAX_ALLOWED_WINDOW_LEN = 5000;
  static constexpr u32 MAX_ALLOWED_REGION_PAD = 1000;

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

  [[nodiscard]] auto BuildWindows() const -> std::vector<WindowPtr>;

 private:
  Params mParams;
  std::unique_ptr<hts::Reference> mRefPtr;

  using ParseRegionResult = hts::Reference::ParseRegionResult;
  absl::flat_hash_set<ParseRegionResult> mInputRegions;

  // Add `regionPadding` to start and end, while checking for coordinate under/over-flow
  void PadInputRegion(ParseRegionResult& result) const;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_WINDOW_BUILDER_H_
