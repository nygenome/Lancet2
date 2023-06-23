#ifndef SRC_LANCET_CLI_CLI_PARAMS_H_
#define SRC_LANCET_CLI_CLI_PARAMS_H_

#include <filesystem>
#include <string>
#include <vector>

#include "lancet/base/types.h"
#include "lancet/caller/variant_call.h"
#include "lancet/core/async_worker.h"
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
  std::filesystem::path mRunStats;
  std::vector<std::string> mInRegions;
  std::vector<std::filesystem::path> mNormalPaths;
  std::vector<std::filesystem::path> mTumorPaths;

  usize mNumWorkerThreads = 2;
  bool mEnableVerboseLogging = false;
#ifndef LANCET_DEVELOP_MODE
  bool mEnableCpuProfiling = false;
#endif

  core::WindowBuilder::Params mWindowBuilder;
  core::VariantBuilder::Params mVariantBuilder;
  // NOLINTEND(misc-non-private-member-variables-in-classes)
};

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_CLI_PARAMS_H_
