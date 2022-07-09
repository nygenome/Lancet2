#include "lancet2/cli.h"

#cmakedefine LANCET2_SANITIZER
#ifndef LANCET2_SANITIZER
#include "mimalloc-new-delete.h"
#endif

auto main(int argc, char** argv) noexcept -> int {
  return lancet2::RunCli(argc, argv);
}
