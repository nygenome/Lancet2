#ifndef SRC_LANCET_CLI_ETA_TIMER_H_
#define SRC_LANCET_CLI_ETA_TIMER_H_

#include "lancet/base/compute_stats.h"
#include "lancet/base/timer.h"
#include "lancet/base/types.h"

namespace lancet::cli {

class EtaTimer {
 public:
  explicit EtaTimer(usize num_iterations);

  void Increment();
  [[nodiscard]] auto EstimatedEta() const -> absl::Duration;
  [[nodiscard]] auto RatePerSecond() const -> f64;

 private:
  usize mNumDone = 0;
  usize mNumTotal = 0;
  Timer mProgressTimer;
  OnlineStats mRunStats;
};

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_ETA_TIMER_H_
