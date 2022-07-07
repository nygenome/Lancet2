#pragma once

#include <string>
#include <string_view>

#include "lancet2/sized_ints.h"

namespace lancet2 {
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-const-variable"
#endif
static constexpr char ALIGN_GAP = '-';
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

struct AlnSeqs {
  std::string ref;
  std::string qry;
};

struct AlnSeqsView {
  std::string_view ref{};
  std::string_view qry{};
};

[[nodiscard]] auto Align(std::string_view ref, std::string_view qry) -> AlnSeqs;

[[nodiscard]] auto TrimEndGaps(AlnSeqsView* aln) -> usize;

[[nodiscard]] auto EdlibAlign(std::string_view ref, std::string_view qry) -> AlnSeqs;
}  // namespace lancet2
