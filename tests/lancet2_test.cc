#define CATCH_CONFIG_RUNNER
#include "catch2/catch_session.hpp"

auto main(int argc, char** argv) -> int {
  // Do any global test setup here
  int result = Catch::Session().run(argc, argv);
  // Do any global test teardown here
  return result;
}
