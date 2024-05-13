#include <iostream>

#include "absl/cleanup/cleanup.h"
#include "absl/debugging/failure_signal_handler.h"
#include "absl/debugging/symbolize.h"
#include "lancet/base/logging.h"
#include "lancet/cli/cli_interface.h"
#include "mimalloc-override.h"  // NOLINT(misc-include-cleaner)
#include "mimalloc.h"
#include "spdlog/spdlog.h"

auto main(const int argc, const char** argv) -> int {
  // Initialize the symbolizer to get a human-readable stack trace
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  absl::InitializeSymbolizer(argv[0]);
  const absl::FailureSignalHandlerOptions options{
      .symbolize_stacktrace = true, .use_alternate_stack = true, .call_previous_handler = false};
  absl::InstallFailureSignalHandler(options);

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
