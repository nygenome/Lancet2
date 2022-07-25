#pragma once

#include <array>
#include <utility>

namespace lancet2 {
struct HpCov {
  HpCov() = default;

  u16 FwdCov = 0;  // NOLINT
  u16 RevCov = 0;  // NOLINT
  u16 HP0 = 0;     // NOLINT
  u16 HP1 = 0;     // NOLINT
  u16 HP2 = 0;     // NOLINT

  [[nodiscard]] auto GetTotalCov() const -> u16 { return FwdCov + RevCov; }
};

struct VariantHpCov {
  VariantHpCov(HpCov ref, HpCov alt) : RefAllele(ref), AltAllele(alt) {}
  VariantHpCov() = default;

  [[nodiscard]] auto TotalRefCov() const -> u16 { return RefAllele.GetTotalCov(); }
  [[nodiscard]] auto TotalAltCov() const -> u16 { return AltAllele.GetTotalCov(); }
  [[nodiscard]] auto TotalCov() const -> u16 { return RefAllele.GetTotalCov() + AltAllele.GetTotalCov(); }

  [[nodiscard]] auto RefHP(Haplotype hp) const -> u16 {
    switch (hp) {
      case Haplotype::UNASSIGNED:
        return RefAllele.HP0;
      case Haplotype::FIRST:
        return RefAllele.HP1;
      case Haplotype::SECOND:
        return RefAllele.HP2;
    }
    __builtin_unreachable();
  }

  [[nodiscard]] auto AltHP(Haplotype hp) const -> u16 {
    switch (hp) {
      case Haplotype::UNASSIGNED:
        return AltAllele.HP0;
      case Haplotype::FIRST:
        return AltAllele.HP1;
      case Haplotype::SECOND:
        return AltAllele.HP2;
    }
    __builtin_unreachable();
  }

  [[nodiscard]] auto TotalHP(Haplotype hp) const -> u16 { return RefHP(hp) + AltHP(hp); }

  HpCov RefAllele;  // NOLINT
  HpCov AltAllele;  // NOLINT

  [[nodiscard]] auto VAF() const -> double {
    const auto altTotal = AltAllele.GetTotalCov();
    if (altTotal == 0) return 0.0;
    return static_cast<double>(altTotal) / static_cast<double>(RefAllele.GetTotalCov() + altTotal);
  }
};
}  // namespace lancet2
