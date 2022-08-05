#include "lancet2/transcript.h"

#include <algorithm>
#include <cmath>

#include "absl/hash/internal/city.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"
#include "lancet2/kmer.h"
#include "lancet2/utils.h"

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

void Transcript::Finalize() {
  if (kind == TranscriptCode::REF_MATCH) {
    isFinalized = true;
    return;
  }

  // get rid of alignment skip chars `-` in ref and alt alleles
  refAllele.erase(std::remove(refAllele.begin(), refAllele.end(), '-'), refAllele.end());
  altAllele.erase(std::remove(altAllele.begin(), altAllele.end(), '-'), altAllele.end());

  const auto refLen = refAllele.length();
  const auto altLen = altAllele.length();
  switch (kind) {
    case TranscriptCode::SNV:
      varLen = 1;
      break;
    case TranscriptCode::INSERTION:
      varLen = altLen;
      break;
    case TranscriptCode::DELETION:
      varLen = refLen;
      break;
    case TranscriptCode::COMPLEX:
      varLen = refLen == altLen ? altLen : refLen > altLen ? (refLen - altLen) : (altLen - refLen);
      break;
    default:
      break;
  }

  // Add previous base for InDels and complex events
  if (kind != TranscriptCode::SNV) {
    refAllele.insert(refAllele.begin(), 1, prevRefBase);
    altAllele.insert(altAllele.begin(), 1, prevAltBase);
    genomeRefPos--;
    idxs.refStart--;
    idxs.altStart--;
    idxs.refEnd--;
    idxs.altEnd--;
  }

  isFinalized = true;
}

auto Transcript::GetAlleleHashes() const noexcept -> AlleleHashes {
  const auto ref = absl::StrFormat("%s_%d_%s_%d", chromName, genomeRefPos, refAllele, varLen);
  const auto alt = absl::StrFormat("%s_%d_%s_%d", chromName, genomeRefPos, altAllele, varLen);

  return AlleleHashes{
      absl::hash_internal::CityHash64WithSeeds(ref.c_str(), ref.length(), utils::PRIME0, utils::PRIME1),
      absl::hash_internal::CityHash64WithSeeds(alt.c_str(), alt.length(), utils::PRIME0, utils::PRIME1)};
}

void Transcript::BuildHaplotypes(std::string_view refSeq, std::string_view altSeq, usize kmerLen) {
  const auto refAlleleLen = refAllele.length();
  const auto altAlleleLen = altAllele.length();

  const auto refHapLen = std::max(refAlleleLen, kmerLen);
  const auto altHapLen = std::max(altAlleleLen, kmerLen);

  const auto isLongRef = refAlleleLen >= kmerLen;
  const auto isLongAlt = altAlleleLen >= kmerLen;

  hapData.clear();
  hapData.reserve(6);

  {
    // REF ALLELE – HaplotypeCentered
    // ----------x----------
    // --------xxxxx--------
    const auto leftFlank =
        isLongRef ? 0 : static_cast<i64>(std::ceil(static_cast<double>(refHapLen - refAlleleLen) / 2.0));

    auto hapStart0 = static_cast<i64>(idxs.refStart) - leftFlank;
    hapStart0 = hapStart0 < 0 ? 0 : hapStart0;
    hapStart0 = (hapStart0 + refHapLen) >= refSeq.length() ? static_cast<i64>(refSeq.length() - refHapLen) : hapStart0;
    const auto hapSeq = absl::ClippedSubstr(refSeq, hapStart0, refHapLen);
    if (hapSeq.length() == refHapLen) {
      hapData.emplace_back(
          HaplotypeData{Kmer::CanonicalSeqHash(hapSeq), refHapLen, static_cast<usize>(leftFlank), Allele::REF});
    }
  }

  {
    // ALT ALLELE – HaplotypeCentered
    // ----------x----------
    // --------xxxxx--------
    const auto leftFlank =
        isLongAlt ? 0 : static_cast<i64>(std::ceil(static_cast<double>(altHapLen - altAlleleLen) / 2.0));

    auto hapStart0 = static_cast<i64>(idxs.altStart) - leftFlank;
    hapStart0 = hapStart0 < 0 ? 0 : hapStart0;
    hapStart0 = (hapStart0 + altHapLen) >= altSeq.length() ? static_cast<i64>(altSeq.length() - altHapLen) : hapStart0;
    const auto hapSeq = absl::ClippedSubstr(altSeq, hapStart0, altHapLen);
    if (hapSeq.length() == altHapLen) {
      hapData.emplace_back(
          HaplotypeData{Kmer::CanonicalSeqHash(hapSeq), altHapLen, static_cast<usize>(leftFlank), Allele::ALT});
    }
  }
}
}  // namespace lancet2
