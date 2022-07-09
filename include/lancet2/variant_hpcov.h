#pragma once

#include "lancet2/base_hpcov.h"
#include "lancet2/fisher_exact.h"

namespace lancet2 {
struct VariantHpCov {
  VariantHpCov(HpCov ref, HpCov alt) : refAl(ref), altAl(alt) {}
  VariantHpCov() = delete;

  [[nodiscard]] auto TotalRefCov() const -> u16 { return refAl.GetTotalCov(); }
  [[nodiscard]] auto TotalAltCov() const -> u16 { return altAl.GetTotalCov(); }
  [[nodiscard]] auto TotalCov() const -> u16 { return refAl.GetTotalCov() + altAl.GetTotalCov(); }

  [[nodiscard]] auto RefHP(Haplotype hp) const -> u16 {
    switch (hp) {
      case Haplotype::UNASSIGNED:
        return refAl.HP0;
      case Haplotype::FIRST:
        return refAl.HP1;
      case Haplotype::SECOND:
        return refAl.HP2;
    }
    __builtin_unreachable();
  }

  [[nodiscard]] auto AltHP(Haplotype hp) const -> u16 {
    switch (hp) {
      case Haplotype::UNASSIGNED:
        return altAl.HP0;
      case Haplotype::FIRST:
        return altAl.HP1;
      case Haplotype::SECOND:
        return altAl.HP2;
    }
    __builtin_unreachable();
  }

  [[nodiscard]] auto TotalHP(Haplotype hp) const -> u16 { return RefHP(hp) + AltHP(hp); }

  HpCov refAl;  // NOLINT
  HpCov altAl;  // NOLINT

  [[nodiscard]] auto VAF() const -> double {
    const auto altTotal = altAl.GetTotalCov();
    if (altTotal == 0) return 0.0;
    return static_cast<double>(altTotal) / static_cast<double>(refAl.GetTotalCov() + altTotal);
  }
};
}  // namespace lancet2
