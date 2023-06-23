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
  auto& [done_counts, waiting_counts] = mCounters;

  moodycamel::ConsumerToken consumer_token(*mInputPtr);
  const moodycamel::ProducerToken producer_token(*mOutputPtr);

  while (waiting_counts->load(std::memory_order_acquire) != 0) {
    // Check if stop is requested for this thread by the RunMain/caller thread
    if (stop_token.stop_requested()) {
      LOG_DEBUG("Quitting AsyncWorker thread {:#x} after processing {} windows", tid, num_done)
      return;
    }

    // Get the next available unprocessed window from the RunMain/caller thread
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (!mInputPtr->try_dequeue(consumer_token, window_ptr)) continue;

    waiting_counts->fetch_add(-1, std::memory_order_release);
    timer.Reset();

    auto variants = mVariantBuilderPtr->ProcessWindow(std::const_pointer_cast<const Window>(window_ptr));
    mVariantStorePtr->AddVariants(std::move(variants));

    const auto status_code = mVariantBuilderPtr->CurrentStatus();
    auto result = Result{.mGenomeIdx = window_ptr->GenomeIndex(), .mRuntime = timer.Runtime(), .mStatus = status_code};
    
    mOutputPtr->enqueue(producer_token, std::move(result));
    done_counts->fetch_add(-1, std::memory_order_release);
    num_done++;
  }

  LOG_DEBUG("Quitting AsyncWorker thread {:#x} after processing {} windows", tid, num_done)
}

}  // namespace lancet::core
