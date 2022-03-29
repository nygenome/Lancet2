#pragma once

#include "absl/random/random.h"
#include "absl/random/seed_sequences.h"
#include "absl/random/uniform_real_distribution.h"

namespace lancet {
class FractionalSampler {
 public:
  explicit FractionalSampler(double needed_fraction) : fractionToKeep(needed_fraction) {}

  [[nodiscard]] auto ShouldSample() -> bool { return (fractionToKeep == 1.0) || distGen(ubrg) <= fractionToKeep; }

 private:
  double fractionToKeep;
  absl::BitGen ubrg{absl::MakeSeedSeq()};
  absl::uniform_real_distribution<double> distGen{0.0, 1.0};
};
}  // namespace lancet
