#include "lancet/statusor.h"

#include <cstdlib>
#include <iostream>

#include "absl/strings/str_format.h"

namespace lancet::internal_statusor {
void Helper::HandleInvalidStatusCtorArg(absl::Status* status) {
  const char* kMessage = "An OK status is not a valid constructor argument to StatusOr<T>";
  std::cerr << absl::StreamFormat("%s\n", kMessage);
  // Fall back to tensorflow::error::INTERNAL.
  *status = absl::InternalError(kMessage);
}

void Helper::Crash(const absl::Status& status) {
  std::cerr << absl::StreamFormat("Attempting to fetch value instead of handling error: %s\n", status.ToString());
  std::exit(EXIT_FAILURE);
}
}  // namespace lancet::internal_statusor
