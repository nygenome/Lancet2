#include "lancet/transcript.h"

#include <cmath>

#include "absl/strings/str_format.h"
#include "lancet/assert_macro.h"

namespace lancet {
Transcript::Transcript(std::string chrom, std::size_t genome_ref_pos, TranscriptCode k, TranscriptOffsets offs,
                       TranscriptBases bases, std::array<SampleCov, 2> covs, bool somatic_status)
    : chromName(std::move(chrom)),
      genomeRefPos(genome_ref_pos),
      kind(k),
      idxs(offs),
      sampleCovs(covs),
      isSomatic(somatic_status) {
  refSeq.push_back(bases.refBase);
  altSeq.push_back(bases.altBase);
  prevRefBase = bases.prevRefBase;
  prevAltBase = bases.prevAltBase;
}

auto Transcript::HasAltCov() const -> bool {
  const auto isSnv = kind == TranscriptCode::SNV;
  const auto acnf = sampleCovs[0].NonZeroMinimum(Allele::ALT, Strand::FWD, isSnv);
  const auto acnr = sampleCovs[0].NonZeroMinimum(Allele::ALT, Strand::REV, isSnv);
  const auto actf = sampleCovs[1].NonZeroMinimum(Allele::ALT, Strand::FWD, isSnv);
  const auto actr = sampleCovs[1].NonZeroMinimum(Allele::ALT, Strand::REV, isSnv);

  return acnf > 0 || acnr > 0 || actf > 0 || actr > 0;
}

auto Transcript::SetRefEndOffset(std::size_t val) -> Transcript& {
  idxs.refEnd = val;
  return *this;
}

auto Transcript::SetAltEndOffset(std::size_t val) -> Transcript& {
  idxs.altEnd = val;
  return *this;
}

auto Transcript::SetCode(TranscriptCode val) -> Transcript& {
  kind = val;
  return *this;
}

auto Transcript::VariantCov(SampleLabel label) const -> VariantHpCov {
  const auto isSnv = kind == TranscriptCode::SNV;
  using U16 = std::uint16_t;

  if (label == SampleLabel::TUMOR) {
    const auto refFwd = isSomatic ? static_cast<U16>(std::ceil(sampleCovs[1].Mean(Allele::REF, Strand::FWD, false)))
                                  : static_cast<U16>(sampleCovs[1].Minimum(Allele::REF, Strand::FWD, false));
    const auto refRev = isSomatic ? static_cast<U16>(std::ceil(sampleCovs[1].Mean(Allele::REF, Strand::REV, false)))
                                  : static_cast<U16>(sampleCovs[1].Minimum(Allele::REF, Strand::REV, false));

    const auto altFwd = static_cast<U16>(sampleCovs[1].Minimum(Allele::ALT, Strand::FWD, isSnv));
    const auto altRev = static_cast<U16>(sampleCovs[1].Minimum(Allele::ALT, Strand::REV, isSnv));

    const auto refHp0 = isSomatic
                            ? static_cast<U16>(std::ceil(sampleCovs[1].Mean(Allele::REF, Haplotype::UNASSIGNED, false)))
                            : static_cast<U16>(sampleCovs[1].Minimum(Allele::REF, Haplotype::UNASSIGNED, false));
    const auto refHp1 = isSomatic
                            ? static_cast<U16>(std::ceil(sampleCovs[1].Mean(Allele::REF, Haplotype::FIRST, false)))
                            : static_cast<U16>(sampleCovs[1].Minimum(Allele::REF, Haplotype::FIRST, false));
    const auto refHp2 = isSomatic
                            ? static_cast<U16>(std::ceil(sampleCovs[1].Mean(Allele::REF, Haplotype::SECOND, false)))
                            : static_cast<U16>(sampleCovs[1].Minimum(Allele::REF, Haplotype::SECOND, false));

    const auto altHp0 = static_cast<U16>(sampleCovs[1].Minimum(Allele::ALT, Haplotype::UNASSIGNED, isSnv));
    const auto altHp1 = static_cast<U16>(sampleCovs[1].Minimum(Allele::ALT, Haplotype::FIRST, isSnv));
    const auto altHp2 = static_cast<U16>(sampleCovs[1].Minimum(Allele::ALT, Haplotype::SECOND, isSnv));

    return VariantHpCov(HpCov(std::make_pair(refFwd, refRev), {refHp0, refHp1, refHp2}),
                        HpCov(std::make_pair(altFwd, altRev), {altHp0, altHp1, altHp2}));
  }

  const auto refFwd = isSomatic ? static_cast<U16>(std::ceil(sampleCovs[0].Mean(Allele::REF, Strand::FWD, false)))
                                : static_cast<U16>(sampleCovs[0].Minimum(Allele::REF, Strand::FWD, false));
  const auto refRev = isSomatic ? static_cast<U16>(std::ceil(sampleCovs[0].Mean(Allele::REF, Strand::REV, false)))
                                : static_cast<U16>(sampleCovs[0].Minimum(Allele::REF, Strand::REV, false));

  const U16 altFwd = static_cast<U16>(sampleCovs[0].Minimum(Allele::ALT, Strand::FWD, isSnv));
  const U16 altRev = static_cast<U16>(sampleCovs[0].Minimum(Allele::ALT, Strand::REV, isSnv));

  const auto refHp0 = isSomatic
                          ? static_cast<U16>(std::ceil(sampleCovs[0].Mean(Allele::REF, Haplotype::UNASSIGNED, false)))
                          : static_cast<U16>(sampleCovs[0].Minimum(Allele::REF, Haplotype::UNASSIGNED, false));
  const auto refHp1 = isSomatic ? static_cast<U16>(std::ceil(sampleCovs[0].Mean(Allele::REF, Haplotype::FIRST, false)))
                                : static_cast<U16>(sampleCovs[0].Minimum(Allele::REF, Haplotype::FIRST, false));
  const auto refHp2 = isSomatic ? static_cast<U16>(std::ceil(sampleCovs[0].Mean(Allele::REF, Haplotype::SECOND, false)))
                                : static_cast<U16>(sampleCovs[0].Minimum(Allele::REF, Haplotype::SECOND, false));

  const U16 altHp0 = static_cast<U16>(sampleCovs[0].Minimum(Allele::ALT, Haplotype::UNASSIGNED, isSnv));
  const U16 altHp1 = static_cast<U16>(sampleCovs[0].Minimum(Allele::ALT, Haplotype::FIRST, isSnv));
  const U16 altHp2 = static_cast<U16>(sampleCovs[0].Minimum(Allele::ALT, Haplotype::SECOND, isSnv));

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
}  // namespace lancet
