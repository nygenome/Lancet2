#include "lancet/hts/uri_utils.h"

extern "C" {
#include "htslib/hfile.h"
}

#include "absl/strings/match.h"
#include "spdlog/fmt/bundled/format.h"

#include <string>
#include <string_view>

namespace lancet::hts {

auto IsCloudUri(std::string_view uri) -> bool {
  return absl::StartsWith(uri, "gs://") ||
         absl::StartsWith(uri, "s3://") ||
         absl::StartsWith(uri, "http://") ||
         absl::StartsWith(uri, "https://") ||
         absl::StartsWith(uri, "ftp://") ||
         absl::StartsWith(uri, "ftps://");
}

auto ValidateCloudAccess(std::string const& uri, std::string const& mode) -> std::string {
  // htslib C FFI requires vararg-style API
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg)
  auto* fptr = hopen(uri.c_str(), mode.c_str());
  if (fptr == nullptr) {
    return fmt::format("Could not open cloud/web resource: {}", uri);
  }
  if (hclose(fptr) < 0) {
    return fmt::format("Failed to close cloud/web resource connection: {}", uri);
  }
  return "";
}

}  // namespace lancet::hts
