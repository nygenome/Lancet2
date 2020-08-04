#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

#include "lancet/base_hpcov.h"
#include "lancet/core_enums.h"
#include "lancet/cov_stats.h"

namespace lancet {
class SampleCov {
 public:
  SampleCov(const BaseHpCov& ref, const BaseHpCov& alt);
  SampleCov() = delete;

  void PushRefAlt(const BaseHpCov& ref, const BaseHpCov& alt);
  void PushRef(const BaseHpCov& ref);
  void PushAlt(const BaseHpCov& alt);

  [[nodiscard]] auto Mean(Allele al, Strand st, bool bqpass) const -> float;
  [[nodiscard]] auto Mean(Allele al, Haplotype hp, bool bqpass) const -> float;

  [[nodiscard]] auto NonZeroMean(Allele al, Strand st, bool bqpass) const -> float;
  [[nodiscard]] auto NonZeroMean(Allele al, Haplotype hp, bool bqpass) const -> float;

  [[nodiscard]] auto Minimum(Allele al, Strand st, bool bqpass) const -> std::uint16_t;
  [[nodiscard]] auto Minimum(Allele al, Haplotype hp, bool bqpass) const -> std::uint16_t;

  [[nodiscard]] auto NonZeroMinimum(Allele al, Strand st, bool bqpass) const -> std::uint16_t;
  [[nodiscard]] auto NonZeroMinimum(Allele al, Haplotype hp, bool bqpass) const -> std::uint16_t;

 private:
  static constexpr std::size_t numStats = 20;
  std::array<CovStats, numStats> data;

  void PushAllele(const BaseHpCov& d, Allele al);

  [[nodiscard]] static auto ToIdx(Allele al, Strand st, bool bqpass) -> std::size_t;
  [[nodiscard]] static auto ToIdx(Allele al, Haplotype hp, bool bqpass) -> std::size_t;

  static constexpr std::size_t refFwdRawPos = 0;
  static constexpr std::size_t refFwdBqPassPos = 1;
  static constexpr std::size_t refRevRawPos = 2;
  static constexpr std::size_t refRevBqPassPos = 3;

  static constexpr std::size_t refRawHp0Pos = 4;
  static constexpr std::size_t refRawHp1Pos = 5;
  static constexpr std::size_t refRawHp2Pos = 6;
  static constexpr std::size_t refBqPassHp0Pos = 7;
  static constexpr std::size_t refBqPassHp1Pos = 8;
  static constexpr std::size_t refBqPassHp2Pos = 9;

  static constexpr std::size_t altFwdRawPos = 10;
  static constexpr std::size_t altFwdBqPassPos = 11;
  static constexpr std::size_t altRevRawPos = 12;
  static constexpr std::size_t altRevBqPassPos = 13;

  static constexpr std::size_t altRawHp0Pos = 14;
  static constexpr std::size_t altRawHp1Pos = 15;
  static constexpr std::size_t altRawHp2Pos = 16;
  static constexpr std::size_t altBqPassHp0Pos = 17;
  static constexpr std::size_t altBqPassHp1Pos = 18;
  static constexpr std::size_t altBqPassHp2Pos = 19;
};
}  // namespace lancet
