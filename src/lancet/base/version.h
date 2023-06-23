#ifndef SRC_LANCET_BASE_VERSION_H_
#define SRC_LANCET_BASE_VERSION_H_

#include "lancet_version.h"
#include "spdlog/fmt/fmt.h"

static constexpr auto LANCET_VERSION_TAG = lancet::VersionTag;
static constexpr auto LANCET_GIT_BRANCH = lancet::GitBranch;
static constexpr auto LANCET_GIT_REVISION = lancet::GitRevision;

[[nodiscard]] inline auto LancetFullVersion() -> std::string {
  static const auto result = fmt::format("{}-{}-{}", LANCET_VERSION_TAG, LANCET_GIT_BRANCH, LANCET_GIT_REVISION);
  return result;
}

#endif  // SRC_LANCET_BASE_VERSION_H_
