#ifndef SRC_LANCET_CLI_PIPELINE_RUNNER_H_
#define SRC_LANCET_CLI_PIPELINE_RUNNER_H_

#include <memory>
#include <string>

#include "lancet/cli/cli_params.h"

namespace lancet::cli {

class PipelineRunner {
 public:
  explicit PipelineRunner(std::shared_ptr<CliParams> params);

  [[noreturn]] void Run();

 private:
  std::shared_ptr<CliParams> mParamsPtr;

  [[nodiscard]] static auto BuildWindows(const CliParams& params) -> std::vector<core::WindowPtr>;
  [[nodiscard]] static auto BuildVcfHeader(const CliParams& params) -> std::string;

  void ValidateAndPopulateParams();
};

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_PIPELINE_RUNNER_H_
