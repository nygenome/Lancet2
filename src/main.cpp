#include "lancet/cli.h"
#include "mimalloc-new-delete.h"

auto main(int argc, char** argv) noexcept -> int {
  mi_option_set(mi_option_show_stats, 1);
  return lancet::RunCli(argc, argv);
}
