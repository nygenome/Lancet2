#include "lancet/core/async_worker.h"

#include "lancet/base/logging.h"
#include "lancet/base/timer.h"
#include "lancet/base/types.h"
#include "lancet/core/window.h"

#include <chrono>
#include <concurrentqueue.h>
#include <memory>
#include <stop_token>
#include <thread>
#include <utility>

namespace lancet::core {

// NOLINTNEXTLINE(performance-unnecessary-value-param)
void AsyncWorker::Process(std::stop_token stop_token,
                          moodycamel::ProducerToken const& /*in_token*/) {
  static thread_local auto const THREAD_ID =
      absl::Hash<std::thread::id>()(std::this_thread::get_id());
  LOG_DEBUG("Starting AsyncWorker thread {:#x}", THREAD_ID)

  Timer timer;
  usize num_done = 0;
  auto window_ptr = std::make_shared<Window>();
  moodycamel::ProducerToken const out_token(*mOutPtr);
  constexpr auto QUEUE_TIMEOUT = std::chrono::milliseconds(10);

  while (true) {
    // Check if stop is requested for this thread by the RunMain/caller thread
    if (stop_token.stop_requested())
      break;

    // Get the next available unprocessed window from the RunMain/caller thread.
    // NOTE: wait_dequeue_timed serves as a blocking futex call preventing hardware thread-spin.
    // If the timeout triggers without a payload, we seamlessly `continue` the loop,
    // structurally allowing the top-level stop_token evaluation to re-poll state.
    if (!mInPtr->wait_dequeue_timed(window_ptr, QUEUE_TIMEOUT))
      continue;

    timer.Reset();
    auto variants = mBuilderPtr->ProcessWindow(std::const_pointer_cast<Window const>(window_ptr));
    mStorePtr->AddVariants(std::move(variants));

    auto const status_code = mBuilderPtr->CurrentStatus();
    mOutPtr->enqueue(out_token, Result{.mGenomeIdx = window_ptr->GenomeIndex(),
                                       .mRuntime = timer.Runtime(),
                                       .mStatus = status_code});
    num_done++;
  }

  LOG_DEBUG("Quitting AsyncWorker thread {:#x} after processing {} windows", THREAD_ID, num_done)
}

}  // namespace lancet::core
