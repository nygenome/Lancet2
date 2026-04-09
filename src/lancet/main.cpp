#include "lancet/base/logging.h"
#include "lancet/cli/cli_interface.h"

#include "absl/cleanup/cleanup.h"

#include <iostream>
#ifndef __APPLE__
#include "mimalloc-override.h"
#endif
#include "mimalloc.h"
#include "spdlog/spdlog.h"

auto main(int const argc, char const** argv) -> int {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  lancet::RegisterLancetLogger();

  absl::Cleanup const cleanup = []() -> void {
    spdlog::shutdown();
    mi_collect(true);
  };

  lancet::cli::CliInterface cli;
  return cli.RunMain(argc, argv);
}
