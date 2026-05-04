#include "lancet/base/version.h"

#include "catch_amalgamated.hpp"

#include <string>

namespace lancet::base::tests {

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  LancetFullVersion                                                       ║
// ║                                                                          ║
// ║  Composed at compile time from \`VERSION_TAG\`, \`GIT_BRANCH\`, and          ║
// ║  \`GIT_REVISION\` constants populated by CMake (lancet_version.h.inc). The ║
// ║  exact branch and revision are environment-dependent, so we assert       ║
// ║  structural invariants — non-empty, version-tag prefix — rather than     ║
// ║  exact strings.                                                          ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("LancetFullVersion returns a non-empty string", "[lancet][base][LancetFullVersion]") {
  // The version string is composed at compile time and must always be
  // populated (at minimum, with VERSION_TAG). An empty result indicates a
  // build-system regression (lancet_version.h.inc not generated).
  CHECK_FALSE(LancetFullVersion().empty());
}

TEST_CASE("LancetFullVersion always begins with the version tag",
          "[lancet][base][LancetFullVersion]") {
  // The implementation either returns VERSION_TAG verbatim (when branch +
  // revision are empty) or "<VERSION_TAG>_<branch>_<revision>". Either way,
  // the string starts with the version tag.
  auto const full = LancetFullVersion();
  std::string const tag = VERSION_TAG;
  REQUIRE(!tag.empty());
  CHECK(full.compare(0, tag.size(), tag) == 0);
}

TEST_CASE("LancetFullVersion preserves git branch and revision when both present",
          "[lancet][base][LancetFullVersion]") {
  // When the build environment populates GIT_BRANCH and GIT_REVISION (the
  // common case for any developer build), the assembled string includes
  // both. We test conditionally: in a stripped-tag-only build the
  // contains-branch claim wouldn't apply.
  std::string const branch = GIT_BRANCH;
  std::string const revision = GIT_REVISION;
  auto const full = LancetFullVersion();

  if (!branch.empty() && !revision.empty()) {
    INFO("full=\"" << full << "\" branch=\"" << branch << "\" rev=\"" << revision << "\"");
    CHECK(full.find(branch) != std::string::npos);
    CHECK(full.find(revision) != std::string::npos);
  }
}

}  // namespace lancet::base::tests
