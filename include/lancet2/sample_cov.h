#pragma once

#include <array>

#include "lancet2/base_hpcov.h"
#include "lancet2/core_enums.h"
#include "lancet2/cov_stats.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class SampleCov {
 public:
  SampleCov(const BaseHpCov& ref, const BaseHpCov& alt);
  SampleCov() = delete;

  void PushRefAlt(const BaseHpCov& ref, const BaseHpCov& alt);
  void PushRef(const BaseHpCov& ref);
  void PushAlt(const BaseHpCov& alt);

  [[nodiscard]] auto GetMean(Allele al, Strand st, bool bqpass) const -> float;
  [[nodiscard]] auto GetMean(Allele al, Haplotype hp, bool bqpass) const -> float;

  [[nodiscard]] auto GetNonZeroMean(Allele al, Strand st, bool bqpass) const -> float;
  [[nodiscard]] auto GetNonZeroMean(Allele al, Haplotype hp, bool bqpass) const -> float;

  [[nodiscard]] auto GetMinimum(Allele al, Strand st, bool bqpass) const -> u16;
  [[nodiscard]] auto GetMinimum(Allele al, Haplotype hp, bool bqpass) const -> u16;

  [[nodiscard]] auto GetNonZeroMinimum(Allele al, Strand st, bool bqpass) const -> u16;
  [[nodiscard]] auto GetNonZeroMinimum(Allele al, Haplotype hp, bool bqpass) const -> u16;

 private:
  static constexpr usize numStats = 20;
  std::array<CovStats, numStats> data;

  void PushAllele(const BaseHpCov& d, Allele al);

  [[nodiscard]] static auto ToIdx(Allele al, Strand st, bool bqpass) -> usize;
  [[nodiscard]] static auto ToIdx(Allele al, Haplotype hp, bool bqpass) -> usize;

  static constexpr usize refFwdRawPos = 0;
  static constexpr usize refFwdBqPassPos = 1;
  static constexpr usize refRevRawPos = 2;
  static constexpr usize refRevBqPassPos = 3;

  static constexpr usize refRawHp0Pos = 4;
  static constexpr usize refRawHp1Pos = 5;
  static constexpr usize refRawHp2Pos = 6;
  static constexpr usize refBqPassHp0Pos = 7;
  static constexpr usize refBqPassHp1Pos = 8;
  static constexpr usize refBqPassHp2Pos = 9;

  static constexpr usize altFwdRawPos = 10;
  static constexpr usize altFwdBqPassPos = 11;
  static constexpr usize altRevRawPos = 12;
  static constexpr usize altRevBqPassPos = 13;

  static constexpr usize altRawHp0Pos = 14;
  static constexpr usize altRawHp1Pos = 15;
  static constexpr usize altRawHp2Pos = 16;
  static constexpr usize altBqPassHp0Pos = 17;
  static constexpr usize altBqPassHp1Pos = 18;
  static constexpr usize altBqPassHp2Pos = 19;
};
}  // namespace lancet2
