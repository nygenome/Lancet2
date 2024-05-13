#ifndef SRC_LANCET_CLI_CLI_PARAMS_H_
#define SRC_LANCET_CLI_CLI_PARAMS_H_

#include <filesystem>
#include <string>
#include <vector>

#include "lancet/base/types.h"
#include "lancet/core/variant_builder.h"
#include "lancet/core/window_builder.h"

namespace lancet::cli {

class CliParams {
 public:
  CliParams() = default;

  // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
  std::string mFullCmdLine;
  std::filesystem::path mOutVcfGz;
  std::filesystem::path mBedFile;
  std::vector<std::string> mInRegions;

  usize mNumWorkerThreads = 2;
  bool mEnableVerboseLogging = false;

  core::WindowBuilder::Params mWindowBuilder;
  core::VariantBuilder::Params mVariantBuilder;
  // NOLINTEND(misc-non-private-member-variables-in-classes)
};

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_CLI_PARAMS_H_
