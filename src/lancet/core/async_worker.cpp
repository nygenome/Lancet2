#include "lancet/core/async_worker.h"

#include <stop_token>

#include "concurrentqueue.h"
#include "lancet/base/logging.h"
#include "lancet/base/timer.h"
#include "lancet/base/types.h"
#include "lancet/core/window.h"

namespace lancet::core {

// NOLINTNEXTLINE(performance-unnecessary-value-param)
void AsyncWorker::Process(std::stop_token stop_token, const moodycamel::ProducerToken& in_token) {
  static thread_local const auto tid = absl::Hash<std::thread::id>()(std::this_thread::get_id());
  LOG_DEBUG("Starting AsyncWorker thread {:#x}", tid)

  Timer timer;
  usize num_done = 0;
  auto window_ptr = std::make_shared<Window>();
  const moodycamel::ProducerToken out_token(*mOutPtr);

  while (true) {
    // Check if stop is requested for this thread by the RunMain/caller thread
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (stop_token.stop_requested()) break;

    // Get the next available unprocessed window from the RunMain/caller thread
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (!mInPtr->try_dequeue_from_producer(in_token, window_ptr)) continue;

    timer.Reset();
    auto variants = mBuilderPtr->ProcessWindow(std::const_pointer_cast<const Window>(window_ptr));
    mStorePtr->AddVariants(std::move(variants));

    const auto status_code = mBuilderPtr->CurrentStatus();
    mOutPtr->enqueue(out_token, Result{window_ptr->GenomeIndex(), timer.Runtime(), status_code});
    num_done++;
  }

  LOG_DEBUG("Quitting AsyncWorker thread {:#x} after processing {} windows", tid, num_done)
}

}  // namespace lancet::core
