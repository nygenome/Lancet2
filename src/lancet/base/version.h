#ifndef SRC_LANCET_BASE_VERSION_H_
#define SRC_LANCET_BASE_VERSION_H_

#include "lancet_version.h"
#include "spdlog/fmt/bundled/format.h"

#include <string>

namespace lancet::base {

static constexpr auto VERSION_TAG = lancet::VersionTag;
static constexpr auto GIT_BRANCH = lancet::GitBranch;
static constexpr auto GIT_REVISION = lancet::GitRevision;

[[nodiscard]] inline auto LancetFullVersion() -> std::string {
  if (GIT_BRANCH[0] == '\0' || GIT_REVISION[0] == '\0') {
    return VERSION_TAG;
  }

  static auto const RESULT = fmt::format("{}_{}_{}", VERSION_TAG, GIT_BRANCH, GIT_REVISION);
  return RESULT;
}

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_VERSION_H_
