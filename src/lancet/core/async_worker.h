#ifndef SRC_LANCET_CORE_ASYNC_WORKER_H_
#define SRC_LANCET_CORE_ASYNC_WORKER_H_

#include <array>
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

  using InQueuePtr = std::shared_ptr<InputQueue>;
  using OutQueuePtr = std::shared_ptr<OutputQueue>;
  using VariantStorePtr = std::shared_ptr<VariantStore>;
  using VariantBuilderPtr = std::unique_ptr<VariantBuilder>;
  using BuilderParamsPtr = std::shared_ptr<const VariantBuilder::Params>;

  AsyncWorker(InQueuePtr in_queue, OutQueuePtr out_queue, VariantStorePtr vstore, BuilderParamsPtr prms)
      : mInPtr(std::move(in_queue)), mOutPtr(std::move(out_queue)), mStorePtr(std::move(vstore)),
        mBuilderPtr(std::make_unique<VariantBuilder>(std::move(prms))) {}

  void Process(std::stop_token stop_token, const moodycamel::ProducerToken& in_token);

 private:
  InQueuePtr mInPtr;
  OutQueuePtr mOutPtr;
  VariantStorePtr mStorePtr;
  VariantBuilderPtr mBuilderPtr;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_ASYNC_WORKER_H_
