#include "lancet2/transcript.h"

#include "absl/strings/str_format.h"

namespace lancet2 {
Transcript::Transcript(std::string chrom, usize genome_ref_pos, TranscriptCode kind, Transcript::Offsets offsets,
                       Transcript::Bases bases)
    : chromName(std::move(chrom)), genomeRefPos(genome_ref_pos), kind(kind), idxs(offsets),
      prevRefBase(bases.prevRefBase), prevAltBase(bases.prevAltBase) {
  refAllele.push_back(bases.refBase);
  altAllele.push_back(bases.altBase);
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

auto Transcript::STRResult() const noexcept -> std::string {
  return strQry.foundSTR ? absl::StrFormat("%d:%s", strQry.strLength, strQry.strMotif) : "";
}

auto Transcript::AddSTRResult(const TandemRepeatResult& val) -> Transcript& {
  strQry = val;
  return *this;
}

auto Transcript::AddRefBase(const char& b) -> Transcript& {
  refAllele.push_back(b);
  return *this;
}

auto Transcript::AddAltBase(const char& b) -> Transcript& {
  altAllele.push_back(b);
  return *this;
}
}  // namespace lancet2
