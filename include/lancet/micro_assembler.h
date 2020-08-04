#pragma once

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/types/span.h"
#include "lancet/cli_params.h"
#include "lancet/notification.h"
#include "lancet/ref_window.h"
#include "lancet/timer.h"
#include "lancet/variant_store.h"

namespace lancet {
class MicroAssembler {
 public:
  using NotificationPtr = std::shared_ptr<Notification<std::size_t>>;
  explicit MicroAssembler(absl::Span<std::shared_ptr<const RefWindow>> ws, absl::Span<NotificationPtr> notis,
                          std::shared_ptr<const CliParams> p)
      : windows(ws), notifiers(notis), params(std::move(p)) {}

  MicroAssembler() = default;

  [[nodiscard]] auto Process(const std::shared_ptr<VariantStore>& store) const -> absl::Status;

 private:
  absl::Span<std::shared_ptr<const RefWindow>> windows;
  absl::Span<NotificationPtr> notifiers;
  std::shared_ptr<const CliParams> params;

  [[nodiscard]] auto ProcessWindow(const std::shared_ptr<const RefWindow>& w,
                                   const std::shared_ptr<VariantStore>& store) const -> absl::Status;

  [[nodiscard]] auto ShouldSkipWindow(const std::shared_ptr<const RefWindow>& w, Timer* T) const -> bool;
};
}  // namespace lancet
