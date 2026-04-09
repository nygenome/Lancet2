#ifndef SRC_LANCET_CLI_PIPELINE_RUNNER_H_
#define SRC_LANCET_CLI_PIPELINE_RUNNER_H_

#include "lancet/cli/cli_params.h"
#include "lancet/core/window.h"
#include "lancet/core/window_builder.h"

#include <memory>
#include <string>
#include <vector>

namespace lancet::cli {

class PipelineRunner {
 public:
  explicit PipelineRunner(std::shared_ptr<CliParams> params);

  [[noreturn]] void Run();

 private:
  std::shared_ptr<CliParams> mParamsPtr;

  // ---------------------------------------------------------------------------
  // Modularized helpers extracted from the monolithic Run() method
  // ---------------------------------------------------------------------------

  /// Initializes the WindowBuilder from CLI params, populates regions, and sorts them.
  [[nodiscard]] static auto InitWindowBuilder(CliParams const& params) -> core::WindowBuilder;

  /// Builds the full VCF header string from CLI params and reference metadata.
  [[nodiscard]] static auto BuildVcfHeader(CliParams const& params) -> std::string;

  /// Validates BAM/CRAM inputs and populates derived params (e.g. MD tag check).
  void ValidateAndPopulateParams();
};

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_PIPELINE_RUNNER_H_
