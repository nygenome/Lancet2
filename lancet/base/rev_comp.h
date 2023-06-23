#ifndef SRC_LANCET_BASE_REV_COMP_H_
#define SRC_LANCET_BASE_REV_COMP_H_

#include <ranges>
#include <string>
#include <string_view>

#include "lancet/base/types.h"

[[nodiscard]] inline auto RevComp(const char& base) -> char {
  switch (base) {
    case 'A':
    case 'a':
      return 'T';

    case 'C':
    case 'c':
      return 'G';

    case 'G':
    case 'g':
      return 'C';

    case 'T':
    case 't':
      return 'A';

    default:
      return 'N';
  }
}

[[nodiscard]] inline auto RevComp(std::string_view seq) -> std::string {
  std::string rev_comp_seq(seq.size(), 'N');
  usize rc_idx = 0;
  for (const char& itr : std::ranges::reverse_view(seq)) {
    rev_comp_seq[rc_idx] = RevComp(itr);
    ++rc_idx;
  }
  return rev_comp_seq;
}

#endif  // SRC_LANCET_BASE_REV_COMP_H_
