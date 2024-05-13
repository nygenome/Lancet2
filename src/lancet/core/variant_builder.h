#ifndef SRC_LANCET_CORE_VARIANT_BUILDER_H_
#define SRC_LANCET_CORE_VARIANT_BUILDER_H_

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "lancet/base/types.h"
#include "lancet/caller/genotyper.h"
#include "lancet/caller/variant_call.h"
#include "lancet/cbdg/graph.h"
#include "lancet/core/read_collector.h"
#include "lancet/core/window.h"

namespace lancet::core {

class VariantBuilder {
 public:
  static constexpr u32 MIN_PHRED_SCORE = 0;
  static constexpr u32 MAX_PHRED_SCORE = 255;

  struct Params {
    bool mSkipActiveRegion = false;
    std::filesystem::path mOutGraphsDir;

    cbdg::Graph::Params mGraphParams;
    ReadCollector::Params mRdCollParams;
  };

  VariantBuilder(std::shared_ptr<const Params> params);

  enum class StatusCode : u8 {
    UNKNOWN = 0,

    SKIPPED_NONLY_REF_BASES = 1,
    SKIPPED_REF_REPEAT_SEEN = 2,
    SKIPPED_INACTIVE_REGION = 3,
    SKIPPED_NOASM_HAPLOTYPE = 4,
    MISSING_NO_MSA_VARIANTS = 5,
    FOUND_GENOTYPED_VARIANT = 6
  };

  [[nodiscard]] auto CurrentStatus() const noexcept -> StatusCode { return mCurrentCode; }

  using WindowResults = std::vector<std::unique_ptr<caller::VariantCall>>;
  [[nodiscard]] auto ProcessWindow(const std::shared_ptr<const Window>& window) -> WindowResults;

 private:
  cbdg::Graph mDebruijnGraph;
  ReadCollector mReadCollector;
  caller::Genotyper mGenotyper;
  std::shared_ptr<const Params> mParamsPtr;
  StatusCode mCurrentCode = StatusCode::UNKNOWN;

  [[nodiscard]] auto MakeGfaPath(const Window& win, usize comp_id) const -> std::filesystem::path;
};

[[nodiscard]] auto ToString(VariantBuilder::StatusCode status_code) -> std::string;

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_VARIANT_BUILDER_H_
