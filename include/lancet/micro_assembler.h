#pragma once

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "lancet/cli_params.h"
#include "lancet/notification.h"
#include "lancet/ref_window.h"
#include "lancet/timer.h"
#include "lancet/variant_store.h"

namespace lancet {
class MicroAssembler {
 public:
  using NotificationPtr = std::shared_ptr<Notification<std::size_t>>;

  explicit MicroAssembler(std::vector<RefWindow>&& ws, std::vector<NotificationPtr>&& notis,
                          std::shared_ptr<const CliParams> p)
      : windows(std::move(ws)), notifiers(std::move(notis)), params(std::move(p)) {}

  MicroAssembler() = default;

  [[nodiscard]] auto Process(const std::shared_ptr<VariantStore>& store) const -> absl::Status;

 private:
  std::vector<RefWindow> windows;
  std::vector<NotificationPtr> notifiers;
  std::shared_ptr<const CliParams> params;

  [[nodiscard]] auto ProcessWindow(const RefWindow& w, const std::shared_ptr<VariantStore>& store) const
      -> absl::Status;

  [[nodiscard]] auto ShouldSkipWindow(const RefWindow& w, Timer* T) const -> bool;
};
}  // namespace lancet
