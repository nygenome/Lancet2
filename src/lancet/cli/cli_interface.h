#ifndef SRC_LANCET_CLI_CLI_INTERFACE_H_
#define SRC_LANCET_CLI_CLI_INTERFACE_H_

#include "lancet/cli/cli_params.h"

#include "CLI/CLI.hpp"

#include <memory>

namespace lancet::cli {

class CliInterface {
 public:
  CliInterface();

  [[nodiscard]] auto RunMain(int argc, char const** argv) -> int;

 private:
  CLI::App mCliApp;
  std::shared_ptr<CliParams> mParamsPtr;

  static void PipelineSubcmd(CLI::App* app, std::shared_ptr<CliParams>& params);
};

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_CLI_INTERFACE_H_
