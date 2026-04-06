#include <iostream>

#include "absl/cleanup/cleanup.h"
#include "lancet/base/logging.h"
#include "lancet/cli/cli_interface.h"
#ifndef __APPLE__
#include "mimalloc-override.h"  // NOLINT(misc-include-cleaner)
#endif
#include "mimalloc.h"
#include "spdlog/spdlog.h"

auto main(const int argc, const char** argv) -> int {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  using namespace lancet;
  RegisterLancetLogger();

  const absl::Cleanup cleanup = [] {
    spdlog::shutdown();
    mi_collect(true);
  };

  cli::CliInterface cli;
  return cli.RunMain(argc, argv);
}
