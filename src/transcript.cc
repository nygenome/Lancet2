#include "lancet2/transcript.h"

#include <cmath>

#include "absl/strings/str_format.h"
#include "lancet2/assert_macro.h"

namespace lancet2 {
Transcript::Transcript(std::string chrom, usize genome_ref_pos, TranscriptCode k, TranscriptOffsets offs,
                       TranscriptBases bases, std::array<SampleCov, 2> covs, bool somatic_status)
    : chromName(std::move(chrom)), genomeRefPos(genome_ref_pos), kind(k), idxs(offs), sampleCovs(covs),
      prevRefBase(bases.prevRefBase), prevAltBase(bases.prevAltBase), isSomatic(somatic_status) {
  refSeq.push_back(bases.refBase);
  altSeq.push_back(bases.altBase);
}

auto Transcript::HasAltCov() const -> bool {
  const auto isSnv = kind == TranscriptCode::SNV;
  const auto acnf = sampleCovs[0].GetNonZeroMinimum(Allele::ALT, Strand::FWD, isSnv);
  const auto acnr = sampleCovs[0].GetNonZeroMinimum(Allele::ALT, Strand::REV, isSnv);
  const auto actf = sampleCovs[1].GetNonZeroMinimum(Allele::ALT, Strand::FWD, isSnv);
  const auto actr = sampleCovs[1].GetNonZeroMinimum(Allele::ALT, Strand::REV, isSnv);

  return acnf > 0 || acnr > 0 || actf > 0 || actr > 0;
}

auto Transcript::SetRefEndOffset(usize val) -> Transcript& {
  idxs.refEnd = val;
  return *this;
}

auto Transcript::SetAltEndOffset(usize val) -> Transcript& {
  idxs.altEnd = val;
  return *this;
}

auto Transcript::SetCode(TranscriptCode val) -> Transcript& {
  kind = val;
  return *this;
}

auto Transcript::VariantCov(SampleLabel label) const -> VariantHpCov {
  const auto isSNV = kind == TranscriptCode::SNV;

  if (label == SampleLabel::TUMOR) {
    const auto refFwd = isSomatic ? static_cast<u16>(std::round(sampleCovs[1].GetMean(Allele::REF, Strand::FWD, false)))
                                  : sampleCovs[1].GetMinimum(Allele::REF, Strand::FWD, false);
    const auto refRev = isSomatic ? static_cast<u16>(std::round(sampleCovs[1].GetMean(Allele::REF, Strand::FWD, false)))
                                  : sampleCovs[1].GetMinimum(Allele::REF, Strand::FWD, false);

    const auto altFwd = sampleCovs[1].GetMinimum(Allele::ALT, Strand::FWD, isSNV);
    const auto altRev = sampleCovs[1].GetMinimum(Allele::ALT, Strand::REV, isSNV);

    const auto refHp0 =
        isSomatic ? static_cast<u16>(std::round(sampleCovs[1].GetMean(Allele::REF, Haplotype::UNASSIGNED, false)))
                  : sampleCovs[1].GetMinimum(Allele::REF, Haplotype::UNASSIGNED, false);
    const auto refHp1 = isSomatic
                            ? static_cast<u16>(std::round(sampleCovs[1].GetMean(Allele::REF, Haplotype::FIRST, false)))
                            : sampleCovs[1].GetMinimum(Allele::REF, Haplotype::FIRST, false);
    const auto refHp2 = isSomatic
                            ? static_cast<u16>(std::round(sampleCovs[1].GetMean(Allele::REF, Haplotype::SECOND, false)))
                            : sampleCovs[1].GetMinimum(Allele::REF, Haplotype::SECOND, false);

    const auto altHp0 = sampleCovs[1].GetMinimum(Allele::ALT, Haplotype::UNASSIGNED, isSNV);
    const auto altHp1 = sampleCovs[1].GetMinimum(Allele::ALT, Haplotype::FIRST, isSNV);
    const auto altHp2 = sampleCovs[1].GetMinimum(Allele::ALT, Haplotype::SECOND, isSNV);

    return VariantHpCov(HpCov(std::make_pair(refFwd, refRev), {refHp0, refHp1, refHp2}),
                        HpCov(std::make_pair(altFwd, altRev), {altHp0, altHp1, altHp2}));
  }

  const auto refFwd = isSomatic ? static_cast<u16>(std::round(sampleCovs[0].GetMean(Allele::REF, Strand::FWD, false)))
                                : sampleCovs[0].GetMinimum(Allele::REF, Strand::FWD, false);
  const auto refRev = isSomatic ? static_cast<u16>(std::round(sampleCovs[0].GetMean(Allele::REF, Strand::FWD, false)))
                                : sampleCovs[0].GetMinimum(Allele::REF, Strand::FWD, false);

  const auto altFwd = isSNV ? sampleCovs[0].GetMinimum(Allele::ALT, Strand::FWD, true)
                            : sampleCovs[0].GetNonZeroMinimum(Allele::ALT, Strand::FWD, false);
  const auto altRev = isSNV ? sampleCovs[0].GetMinimum(Allele::ALT, Strand::REV, true)
                            : sampleCovs[0].GetNonZeroMinimum(Allele::ALT, Strand::REV, false);

  const auto refHp0 =
      isSomatic ? static_cast<u16>(std::round(sampleCovs[0].GetMean(Allele::REF, Haplotype::UNASSIGNED, false)))
                : sampleCovs[0].GetMinimum(Allele::REF, Haplotype::UNASSIGNED, false);
  const auto refHp1 = isSomatic
                          ? static_cast<u16>(std::round(sampleCovs[0].GetMean(Allele::REF, Haplotype::FIRST, false)))
                          : sampleCovs[0].GetMinimum(Allele::REF, Haplotype::FIRST, false);
  const auto refHp2 = isSomatic
                          ? static_cast<u16>(std::round(sampleCovs[0].GetMean(Allele::REF, Haplotype::SECOND, false)))
                          : sampleCovs[0].GetMinimum(Allele::REF, Haplotype::SECOND, false);

  const auto altHp0 = sampleCovs[0].GetMinimum(Allele::ALT, Haplotype::UNASSIGNED, isSNV);
  const auto altHp1 = sampleCovs[0].GetMinimum(Allele::ALT, Haplotype::FIRST, isSNV);
  const auto altHp2 = sampleCovs[0].GetMinimum(Allele::ALT, Haplotype::SECOND, isSNV);

  return VariantHpCov(HpCov(std::make_pair(refFwd, refRev), {refHp0, refHp1, refHp2}),
                      HpCov(std::make_pair(altFwd, altRev), {altHp0, altHp1, altHp2}));
}

auto Transcript::STRResult() const noexcept -> std::string {
  return strQry.foundSTR ? absl::StrFormat("%d:%s", strQry.strLength, strQry.strMotif) : "";
}

auto Transcript::AddCov(SampleLabel label, Allele al, const BaseHpCov& c) -> Transcript& {
  if (label == SampleLabel::TUMOR) {
    al == Allele::REF ? sampleCovs[1].PushRef(c) : sampleCovs[1].PushAlt(c);
    return *this;
  }

  LANCET_ASSERT(label == SampleLabel::NORMAL);  // NOLINT
  al == Allele::REF ? sampleCovs[0].PushRef(c) : sampleCovs[0].PushAlt(c);
  return *this;
}

auto Transcript::AddSTRResult(const TandemRepeatResult& val) -> Transcript& {
  strQry = val;
  return *this;
}

auto Transcript::AddRefBase(const char& b) -> Transcript& {
  refSeq.push_back(b);
  return *this;
}

auto Transcript::AddAltBase(const char& b) -> Transcript& {
  altSeq.push_back(b);
  return *this;
}

auto Transcript::SetSomaticStatus(bool val) -> Transcript& {
  isSomatic = val;
  return *this;
}

auto Transcript::ComputeState() const -> VariantState {
  const auto nmlAlt = VariantCov(SampleLabel::NORMAL).TotalAltCov();
  const auto tmrAlt = VariantCov(SampleLabel::TUMOR).TotalAltCov();

  if (nmlAlt > 0 && tmrAlt > 0) return VariantState::SHARED;
  if (nmlAlt == 0 && tmrAlt > 0) return VariantState::SOMATIC;
  if (nmlAlt > 0 && tmrAlt == 0) return VariantState::NORMAL;
  return VariantState::NONE;
}
}  // namespace lancet2
