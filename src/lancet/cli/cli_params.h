#ifndef SRC_LANCET_CLI_CLI_PARAMS_H_
#define SRC_LANCET_CLI_CLI_PARAMS_H_

#include "lancet/base/types.h"
#include "lancet/core/variant_builder.h"
#include "lancet/core/window_builder.h"

#include <filesystem>
#include <string>
#include <vector>

namespace lancet::cli {

class CliParams {
 public:
  CliParams() = default;

  std::string mFullCmdLine;
  std::filesystem::path mOutVcfGz;
  std::filesystem::path mBedFile;
  std::vector<std::string> mInRegions;

  usize mNumWorkerThreads = 2;
  bool mEnableVerboseLogging = false;

  core::WindowBuilder::Params mWindowBuilder;
  core::VariantBuilder::Params mVariantBuilder;
};

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_CLI_PARAMS_H_
