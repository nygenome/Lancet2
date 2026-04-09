#ifndef SRC_LANCET_BASE_VERSION_H_
#define SRC_LANCET_BASE_VERSION_H_

#include "lancet_version.h"
#include "spdlog/fmt/bundled/core.h"

#include <string>

static constexpr auto LANCET_VERSION_TAG = lancet::VersionTag;
static constexpr auto LANCET_GIT_BRANCH = lancet::GitBranch;
static constexpr auto LANCET_GIT_REVISION = lancet::GitRevision;

[[nodiscard]] inline auto LancetFullVersion() -> std::string {
  static auto const RESULT =
      fmt::format("{}-{}-{}", LANCET_VERSION_TAG, LANCET_GIT_BRANCH, LANCET_GIT_REVISION);
  return RESULT;
}

#endif  // SRC_LANCET_BASE_VERSION_H_
