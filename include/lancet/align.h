#pragma once

#include <cstddef>
#include <string>
#include <string_view>

namespace lancet {
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-const-variable"
#endif
static constexpr char ALIGN_GAP = '-';
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

struct AlignedSequences {
  std::string ref;
  std::string qry;
};

struct AlignedSequencesView {
  std::string_view ref{};
  std::string_view qry{};
};

[[nodiscard]] auto Align(std::string_view ref, std::string_view qry) -> AlignedSequences;

[[nodiscard]] auto TrimEndGaps(AlignedSequencesView* aln) -> std::size_t;
}  // namespace lancet
