#include "lancet2/sample_cov.h"

namespace lancet2 {
SampleCov::SampleCov(const BaseHpCov& ref, const BaseHpCov& alt) { PushRefAlt(ref, alt); }

void SampleCov::PushRefAlt(const BaseHpCov& ref, const BaseHpCov& alt) {
  PushAllele(ref, Allele::REF);
  PushAllele(alt, Allele::ALT);
}

void SampleCov::PushRef(const BaseHpCov& ref) { return PushAllele(ref, Allele::REF); }
void SampleCov::PushAlt(const BaseHpCov& alt) { return PushAllele(alt, Allele::ALT); }

auto SampleCov::GetMean(Allele al, Strand st, bool bqpass) const -> float {
  return data.at(ToIdx(al, st, bqpass)).GetMean();
}
auto SampleCov::GetMean(Allele al, Haplotype hp, bool bqpass) const -> float {
  return data.at(ToIdx(al, hp, bqpass)).GetMean();
}

auto SampleCov::GetNonZeroMean(Allele al, Strand st, bool bqpass) const -> float {
  return data.at(ToIdx(al, st, bqpass)).GetNonZeroMean();
}
auto SampleCov::GetNonZeroMean(Allele al, Haplotype hp, bool bqpass) const -> float {
  return data.at(ToIdx(al, hp, bqpass)).GetNonZeroMean();
}

auto SampleCov::GetMinimum(Allele al, Strand st, bool bqpass) const -> u16 {
  return data.at(ToIdx(al, st, bqpass)).GetMinimum();
}
auto SampleCov::GetMinimum(Allele al, Haplotype hp, bool bqpass) const -> u16 {
  return data.at(ToIdx(al, hp, bqpass)).GetMinimum();
}

auto SampleCov::GetNonZeroMinimum(Allele al, Strand st, bool bqpass) const -> u16 {
  return data.at(ToIdx(al, st, bqpass)).GetNonZeroMinimum();
}
auto SampleCov::GetNonZeroMinimum(Allele al, Haplotype hp, bool bqpass) const -> u16 {
  return data.at(ToIdx(al, hp, bqpass)).GetNonZeroMinimum();
}

void SampleCov::PushAllele(const BaseHpCov& d, Allele al) {
  data.at(ToIdx(al, Strand::FWD, false)).Push(d.raw.fwdCov);
  data.at(ToIdx(al, Strand::REV, false)).Push(d.raw.revCov);

  data.at(ToIdx(al, Strand::FWD, true)).Push(d.bqPass.fwdCov);
  data.at(ToIdx(al, Strand::REV, true)).Push(d.bqPass.revCov);

  data.at(ToIdx(al, Haplotype::UNASSIGNED, false)).Push(d.raw.HP0);
  data.at(ToIdx(al, Haplotype::FIRST, false)).Push(d.raw.HP1);
  data.at(ToIdx(al, Haplotype::SECOND, false)).Push(d.raw.HP2);

  data.at(ToIdx(al, Haplotype::UNASSIGNED, true)).Push(d.bqPass.HP0);
  data.at(ToIdx(al, Haplotype::FIRST, true)).Push(d.bqPass.HP1);
  data.at(ToIdx(al, Haplotype::SECOND, true)).Push(d.bqPass.HP2);
}

auto SampleCov::ToIdx(Allele al, Strand st, bool bqpass) -> usize {
  if (al == Allele::REF) {
    if (st == Strand::FWD) return bqpass ? refFwdBqPassPos : refFwdRawPos;
    return bqpass ? refRevBqPassPos : refRevRawPos;
  }

  if (st == Strand::FWD) return bqpass ? altFwdBqPassPos : altFwdRawPos;
  return bqpass ? altRevBqPassPos : altRevRawPos;
}

auto SampleCov::ToIdx(Allele al, Haplotype hp, bool bqpass) -> usize {
  if (al == Allele::REF) {
    if (bqpass) {
      return hp == Haplotype::UNASSIGNED ? refBqPassHp0Pos : hp == Haplotype::FIRST ? refBqPassHp1Pos : refBqPassHp2Pos;
    }

    return hp == Haplotype::UNASSIGNED ? refRawHp0Pos : hp == Haplotype::FIRST ? refRawHp1Pos : refRawHp2Pos;
  }

  if (bqpass) {
    return hp == Haplotype::UNASSIGNED ? altBqPassHp0Pos : hp == Haplotype::FIRST ? altBqPassHp1Pos : altBqPassHp2Pos;
  }
  return hp == Haplotype::UNASSIGNED ? altRawHp0Pos : hp == Haplotype::FIRST ? altRawHp1Pos : altRawHp2Pos;
}
}  // namespace lancet2
