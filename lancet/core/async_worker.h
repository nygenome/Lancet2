#ifndef SRC_LANCET_CORE_ASYNC_WORKER_H_
#define SRC_LANCET_CORE_ASYNC_WORKER_H_

#include <array>
#include <atomic>
#include <memory>
#include <stop_token>
#include <utility>
#include <vector>

#include "absl/time/time.h"
#include "concurrentqueue.h"
#include "lancet/base/types.h"
#include "lancet/core/variant_builder.h"
#include "lancet/core/variant_store.h"
#include "lancet/core/window.h"

namespace lancet::core {

class AsyncWorker {
 public:
  struct Result {
    // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
    usize mGenomeIdx = 0;
    absl::Duration mRuntime = absl::ZeroDuration();
    VariantBuilder::StatusCode mStatus = VariantBuilder::StatusCode::UNKNOWN;
    // NOLINTEND(misc-non-private-member-variables-in-classes)
  };

  using InputQueue = moodycamel::ConcurrentQueue<WindowPtr>;
  using OutputQueue = moodycamel::ConcurrentQueue<Result>;

  using Input = std::shared_ptr<InputQueue>;
  using Output = std::shared_ptr<OutputQueue>;
  using Store = std::shared_ptr<VariantStore>;
  using Builder = std::unique_ptr<VariantBuilder>;
  using Params = std::shared_ptr<const VariantBuilder::Params>;

  using AtomicCounter = std::shared_ptr<std::atomic_size_t>;
  using DoneAndWaitingCounts = std::array<AtomicCounter, 2>;

  AsyncWorker(Input in_queue, Output out_queue, Store vstore, Params prms, DoneAndWaitingCounts cntrs)
      : mInputPtr(std::move(in_queue)), mOutputPtr(std::move(out_queue)), mVariantStorePtr(std::move(vstore)),
        mVariantBuilderPtr(std::make_unique<VariantBuilder>(std::move(prms))), mCounters(std::move(cntrs)) {}

  void Process(std::stop_token stop_token);

 private:
  Input mInputPtr;
  Output mOutputPtr;
  Store mVariantStorePtr;
  Builder mVariantBuilderPtr;
  DoneAndWaitingCounts mCounters;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_ASYNC_WORKER_H_
