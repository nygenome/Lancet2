#include <memory>

#include "lancet2/cli_params.h"

namespace lancet2 {
[[noreturn]] void RunPipeline(std::shared_ptr<CliParams> params);
}  // namespace lancet2
