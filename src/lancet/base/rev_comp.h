#ifndef SRC_LANCET_BASE_REV_COMP_H_
#define SRC_LANCET_BASE_REV_COMP_H_

#include "lancet/base/types.h"

#include <array>
#include <ranges>
#include <string>
#include <string_view>

// Constexpr lookup table for DNA complement: A↔T, C↔G, else→N.
constexpr auto MakeDnaComplementTable() -> std::array<char, 256> {
  std::array<char, 256> tbl{};
  for (auto& val : tbl) {
    val = 'N';
  }
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
  tbl['A'] = 'T';
  tbl['a'] = 't';
  tbl['T'] = 'A';
  tbl['t'] = 'a';
  tbl['C'] = 'G';
  tbl['c'] = 'g';
  tbl['G'] = 'C';
  tbl['g'] = 'c';
  tbl['N'] = 'N';
  tbl['n'] = 'n';
  // NOLINTEND(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
  return tbl;
}

inline constexpr std::array<char, 256> DNA_COMPLEMENT_TABLE = MakeDnaComplementTable();

[[nodiscard]] inline auto RevComp(char const& base) -> char {
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
  return DNA_COMPLEMENT_TABLE[static_cast<u8>(base)];
}

[[nodiscard]] inline auto RevComp(std::string_view seq) -> std::string {
  std::string rev_comp_seq(seq.size(), 'N');
  usize rc_idx = 0;
  for (char const& itr : std::ranges::reverse_view(seq)) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
    rev_comp_seq[rc_idx] = DNA_COMPLEMENT_TABLE[static_cast<u8>(itr)];
    ++rc_idx;
  }
  return rev_comp_seq;
}

#endif  // SRC_LANCET_BASE_REV_COMP_H_
