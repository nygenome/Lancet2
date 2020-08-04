#include <memory>

#include "lancet/cli_params.h"

namespace lancet {
[[noreturn]] void RunPipeline(std::shared_ptr<CliParams> params);
}  // namespace lancet
