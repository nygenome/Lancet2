#include "lancet/statusor.h"

#include <stdexcept>
#include <string_view>

#include "absl/strings/str_format.h"
#include "spdlog/spdlog.h"

namespace lancet::internal_statusor {
void Helper::HandleInvalidStatusCtorArg(absl::Status* status) {
  static constexpr std::string_view kMessage = "An OK status is not a valid constructor argument to StatusOr<T>";
  SPDLOG_ERROR(kMessage);
  // Fall back to tensorflow::error::INTERNAL.
  *status = absl::InternalError(kMessage);
}

void Helper::Crash(const absl::Status& status) {
  const auto msg = absl::StrFormat("Attempting to fetch value instead of handling error: %s", status.ToString());
  SPDLOG_CRITICAL(msg);
  throw std::runtime_error(msg);
}
}  // namespace lancet::internal_statusor
