#ifndef SRC_LANCET_CLI_CLI_INTERFACE_H_
#define SRC_LANCET_CLI_CLI_INTERFACE_H_

#include <memory>

#include "CLI/CLI.hpp"
#include "lancet/cli/cli_params.h"

namespace lancet::cli {

class CliInterface {
 public:
  CliInterface();

  [[nodiscard]] auto RunMain(int argc, const char** argv) -> int;

 private:
  CLI::App mCliApp;
  std::shared_ptr<CliParams> mParamsPtr;

  static void PipelineSubcmd(CLI::App* app, std::shared_ptr<CliParams>& params);
};

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_CLI_INTERFACE_H_
