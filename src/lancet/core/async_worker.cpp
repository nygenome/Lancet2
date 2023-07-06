#include "lancet/core/async_worker.h"

#include <thread>

#include "absl/hash/hash.h"
#include "lancet/base/logging.h"
#include "lancet/base/repeat.h"
#include "lancet/base/sliding.h"
#include "lancet/base/timer.h"
#include "lancet/caller/msa_builder.h"
#include "lancet/caller/variant_set.h"

namespace lancet::core {

// NOLINTNEXTLINE(performance-unnecessary-value-param)
void AsyncWorker::Process(std::stop_token stop_token) {
  static thread_local const auto tid = absl::Hash<std::thread::id>()(std::this_thread::get_id());
  LOG_DEBUG("Starting AsyncWorker thread {:#x}", tid)

  Timer timer;
  usize num_done = 0;
  auto window_ptr = std::make_shared<Window>();
  moodycamel::ConsumerToken consumer_token(*mInPtr);
  const moodycamel::ProducerToken producer_token(*mOutPtr);

  while (true) {
    // Check if stop is requested for this thread by the RunMain/caller thread
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (stop_token.stop_requested()) break;

    // Get the next available unprocessed window from the RunMain/caller thread
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (!mInPtr->try_dequeue(consumer_token, window_ptr)) continue;

    timer.Reset();
    auto variants = mBuilderPtr->ProcessWindow(std::const_pointer_cast<const Window>(window_ptr));
    mStorePtr->AddVariants(std::move(variants));

    const auto status_code = mBuilderPtr->CurrentStatus();
    mOutPtr->enqueue(producer_token, Result{window_ptr->GenomeIndex(), timer.Runtime(), status_code});
    num_done++;
  }

  LOG_DEBUG("Quitting AsyncWorker thread {:#x} after processing {} windows", tid, num_done)
}

}  // namespace lancet::core
