#pragma once

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/types/span.h"
#include "concurrentqueue.h"
#include "lancet/cli_params.h"
#include "lancet/notification.h"
#include "lancet/ref_window.h"
#include "lancet/timer.h"
#include "lancet/variant_store.h"

namespace lancet {
using WindowPtr = std::shared_ptr<RefWindow>;
using ResultNotifierPtr = std::shared_ptr<Notification<std::size_t>>;
using WindowQueue = moodycamel::ConcurrentQueue<WindowPtr>;

class MicroAssembler {
 public:
  explicit MicroAssembler(std::shared_ptr<WindowQueue> qptr, absl::Span<ResultNotifierPtr> notis,
                          std::shared_ptr<const CliParams> p)
      : windowQPtr(std::move(qptr)), resultNotifiers(notis), params(std::move(p)) {}

  MicroAssembler() = default;

  void Process(const std::shared_ptr<VariantStore>& store) const;

 private:
  std::shared_ptr<WindowQueue> windowQPtr;
  absl::Span<ResultNotifierPtr> resultNotifiers;
  std::shared_ptr<const CliParams> params;

  [[nodiscard]] auto ProcessWindow(const std::shared_ptr<const RefWindow>& w,
                                   const std::shared_ptr<VariantStore>& store) const -> absl::Status;

  [[nodiscard]] auto ShouldSkipWindow(const std::shared_ptr<const RefWindow>& w) const -> bool;
};
}  // namespace lancet
