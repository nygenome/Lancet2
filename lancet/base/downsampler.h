#ifndef SRC_LANCET_BASE_DOWNSAMPLER_H_
#define SRC_LANCET_BASE_DOWNSAMPLER_H_

#include <random>

#include "lancet/base/types.h"

class Downsampler {
 public:
  Downsampler() = default;
  explicit Downsampler(const f64 percent_needed) : mPercentToKeep(percent_needed) {}

  void SetPercentToSample(const f64 percent_needed) {
    mPercentToKeep = percent_needed;
    mPctGenerator.reset();
  }

  [[nodiscard]] auto ShouldSample() -> bool {
    return (mPercentToKeep == 100.0) || mPctGenerator(mRandEngine) <= mPercentToKeep;
  }

 private:
  static constexpr u64 DEFAULT_SEED = 0xa0761d6478bd642fULL;

  f64 mPercentToKeep = 100.0;
  // NOLINTNEXTLINE(cert-msc32-c,cert-msc51-cpp)
  std::mt19937_64 mRandEngine{DEFAULT_SEED};
  std::uniform_real_distribution<f64> mPctGenerator{0.0, 100.0};
};

#endif  // SRC_LANCET_BASE_DOWNSAMPLER_H_
