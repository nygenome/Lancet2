#ifndef SRC_LANCET_CLI_REMAINING_TIMER_H_
#define SRC_LANCET_CLI_REMAINING_TIMER_H_

#include "absl/time/time.h"
#include "lancet/base/compute_stats.h"
#include "lancet/base/types.h"

namespace lancet::cli {

class RemainingTimer {
 public:
  RemainingTimer(usize num_iterations, f64 prior_mean, f64 prior_sd);

  void Update(const absl::Duration& loop_time);
  [[nodiscard]] auto EstimateRemaining() const -> absl::Duration;

 private:
  usize mNumDone = 0;
  usize mNumTotal = 0;
  OnlineStats mRunStats;
};

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_REMAINING_TIMER_H_
