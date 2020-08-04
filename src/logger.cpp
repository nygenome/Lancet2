#include "lancet/logger.h"

#include "absl/time/clock.h"
#include "absl/time/time.h"  // NOLINT

namespace lancet {
auto RFC3339Time() -> std::string {
  return absl::FormatTime(absl::RFC3339_sec, absl::Now(), absl::LocalTimeZone());  // NOLINT
}
}  // namespace lancet
