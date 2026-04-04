#ifndef SRC_LANCET_BASE_REV_COMP_H_
#define SRC_LANCET_BASE_REV_COMP_H_

#include <array>
#include <ranges>
#include <string>
#include <string_view>

#include "lancet/base/types.h"

// Constexpr lookup table for DNA complement: A↔T, C↔G, else→N.
inline constexpr auto MakeDnaComplementTable() -> std::array<char, 256> {
  std::array<char, 256> tbl{};
  for (auto& v : tbl) v = 'N';
  tbl['A'] = 'T';  tbl['a'] = 't';
  tbl['T'] = 'A';  tbl['t'] = 'a';
  tbl['C'] = 'G';  tbl['c'] = 'g';
  tbl['G'] = 'C';  tbl['g'] = 'c';
  tbl['N'] = 'N';  tbl['n'] = 'n';
  return tbl;
}

inline constexpr std::array<char, 256> kDnaComplementTable = MakeDnaComplementTable();

[[nodiscard]] inline auto RevComp(const char& base) -> char {
  return kDnaComplementTable[static_cast<u8>(base)];
}

[[nodiscard]] inline auto RevComp(std::string_view seq) -> std::string {
  std::string rev_comp_seq(seq.size(), 'N');
  usize rc_idx = 0;
  for (const char& itr : std::ranges::reverse_view(seq)) {
    rev_comp_seq[rc_idx] = kDnaComplementTable[static_cast<u8>(itr)];
    ++rc_idx;
  }
  return rev_comp_seq;
}

#endif  // SRC_LANCET_BASE_REV_COMP_H_
