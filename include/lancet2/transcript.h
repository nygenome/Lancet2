#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "absl/types/span.h"
#include "lancet2/core_enums.h"
#include "lancet2/kmer.h"
#include "lancet2/sized_ints.h"
#include "lancet2/tandem_repeat.h"
#include "lancet2/variant_hpcov.h"

namespace lancet2 {
class Transcript {
 public:
  /// 0-based offsets within reference and path sequences of the transcript
  struct Offsets {
    usize refStart = 0;
    usize altStart = 0;
    usize refEnd = 0;
    usize altEnd = 0;
  };

  struct Bases {
    char refBase = '!';
    char altBase = '!';
    char prevRefBase = '!';
    char prevAltBase = '!';
  };

  Transcript(std::string chrom, usize genome_ref_pos, TranscriptCode kind, Transcript::Offsets offsets,
             Transcript::Bases bases);
  Transcript() = delete;

  [[nodiscard]] auto ChromName() const noexcept -> std::string { return chromName; }
  [[nodiscard]] auto Position() const noexcept -> usize { return genomeRefPos; }
  [[nodiscard]] auto RefStartOffset() const noexcept -> usize { return idxs.refStart; }
  [[nodiscard]] auto AltStartOffset() const noexcept -> usize { return idxs.altStart; }

  [[nodiscard]] auto RefEndOffset() const noexcept -> usize { return idxs.refEnd; }
  auto SetRefEndOffset(usize val) -> Transcript&;

  [[nodiscard]] auto AltEndOffset() const noexcept -> usize { return idxs.altEnd; }
  auto SetAltEndOffset(usize val) -> Transcript&;

  [[nodiscard]] auto Code() const noexcept -> TranscriptCode { return kind; }
  auto SetCode(TranscriptCode val) -> Transcript&;

  [[nodiscard]] auto STRResult() const noexcept -> std::string;
  auto AddSTRResult(const TandemRepeatResult& val) -> Transcript&;

  [[nodiscard]] auto RefSeq() const noexcept -> std::string { return refAllele; }
  [[nodiscard]] auto GetRefLength() const noexcept -> usize { return refAllele.length(); }
  auto AddRefBase(const char& b) -> Transcript&;

  [[nodiscard]] auto AltSeq() const noexcept -> std::string { return altAllele; }
  [[nodiscard]] auto GetAltLength() const noexcept -> usize { return altAllele.length(); }
  auto AddAltBase(const char& b) -> Transcript&;

  [[nodiscard]] auto IsFinalized() const noexcept -> bool { return isFinalized; }
  [[nodiscard]] auto GetVariantLength() const noexcept -> usize { return varLen; }
  void Finalize();

  struct AlleleHashes {
    u64 RefHash;
    u64 AltHash;
  };
  [[nodiscard]] auto GetAlleleHashes() const noexcept -> AlleleHashes;

  struct HaplotypeData {
    u64 hapHash;
    usize hapLen;
    usize hapLeftFlank;
    Allele alleleKind;
  };
  void BuildHaplotypes(std::string_view refSeq, std::string_view altSeq, usize kmerLen);

  [[nodiscard]] auto GetHaplotypesData() const -> absl::Span<const HaplotypeData> {
    return absl::MakeConstSpan(hapData);
  }

 private:
  std::string chromName;
  usize genomeRefPos = 0;  // 1-based genome position for VCF
  TranscriptCode kind = TranscriptCode::REF_MATCH;
  Transcript::Offsets idxs;
  TandemRepeatResult strQry;
  std::string refAllele;
  std::string altAllele;
  char prevRefBase = 'N';
  char prevAltBase = 'N';

  bool isFinalized = false;
  usize varLen = 0;
  std::vector<HaplotypeData> hapData;
};
}  // namespace lancet2
