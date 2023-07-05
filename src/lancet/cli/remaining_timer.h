#ifndef SRC_LANCET_CLI_REMAINING_TIMER_H_
#define SRC_LANCET_CLI_REMAINING_TIMER_H_

#include "absl/time/time.h"
#include "lancet/base/compute_stats.h"
#include "lancet/base/types.h"

namespace lancet::cli {

class RemainingTimer {
 public:
  explicit RemainingTimer(usize num_iterations, usize num_threads);

  void Update(const absl::Duration& loop_time);
  [[nodiscard]] auto EstimateRemaining() const -> absl::Duration;
  [[nodiscard]] auto MeanRatePerSecond() const -> f64;

 private:
  usize mNumDone = 0;
  usize mNumTotal = 0;
  usize mNumThreads = 1;
  OnlineStats mRunStats;
};

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_REMAINING_TIMER_H_
