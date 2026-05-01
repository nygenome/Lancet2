#include "lancet/base/crash_handler.h"
#include "lancet/base/logging.h"
#include "lancet/cli/cli_interface.h"

#include "absl/cleanup/cleanup.h"
#ifndef LANCET_SANITIZE_BUILD
#include "mimalloc.h"  // IWYU pragma: keep
#endif
#include "spdlog/sinks/ansicolor_sink.h"
#include "spdlog/spdlog.h"

#if !defined(__APPLE__) && !defined(LANCET_SANITIZE_BUILD)
#include "mimalloc-override.h"  // IWYU pragma: keep
#endif

#include <exception>
#include <iostream>
#include <memory>

#include <cstdio>
#include <cstdlib>

auto main(int const argc, char const** argv) -> int {
  try {
#ifndef LANCET_SANITIZE_BUILD
    // Disable mimalloc's arena purge cycle entirely. Without this, mimalloc
    // calls madvise(MADV_DONTNEED) to decommit freed pages every ~100ms,
    // costing ~200s of CPU time over a full run. Since graph node count is
    // bounded, per-thread heap memory plateaus early and is reused across
    // windows — purging just adds decommit/recommit churn with no benefit.
    // mi_collect(true) at exit releases everything.
    mi_option_set(mi_option_purge_delay, -1);
#endif

    lancet::base::InstallCrashHandler();
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    lancet::RegisterLancetLogger();
    absl::Cleanup const cleanup = []() -> void {
      spdlog::shutdown();
#ifndef LANCET_SANITIZE_BUILD
      mi_collect(true);
#endif
    };

    lancet::cli::CliInterface cli;
    return cli.RunMain(argc, argv);
  } catch (std::exception const& ex) {
    // spdlog may already be shut down via Cleanup unwind, so write directly to stderr.
    fmt::print(stderr, "FATAL: Unhandled exception in main: {}\n", ex.what());
    return EXIT_FAILURE;
  } catch (...) {
    fmt::print(stderr, "FATAL: Unhandled exception of unknown type in main\n");
    return EXIT_FAILURE;
  }
}
