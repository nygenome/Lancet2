#ifndef SRC_LANCET_HTS_PHRED_QUALITY_H_
#define SRC_LANCET_HTS_PHRED_QUALITY_H_

#include "lancet/base/types.h"

namespace lancet::hts {

static constexpr u8 MAX_PHRED_SCORE = 255;
[[nodiscard]] auto PhredToErrorProb(u32 phred_score) -> f64;
[[nodiscard]] auto ErrorProbToPhred(f64 prob) -> u8;

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_PHRED_QUALITY_H_
