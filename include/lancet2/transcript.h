#pragma once

#include <array>
#include <string>
#include <string_view>

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
  auto AddRefBase(const char& b) -> Transcript&;

  [[nodiscard]] auto AltSeq() const noexcept -> std::string { return altAllele; }
  auto AddAltBase(const char& b) -> Transcript&;

  [[nodiscard]] auto PrevRefBase() const noexcept -> char { return prevRefBase; }
  [[nodiscard]] auto PrevAltBase() const noexcept -> char { return prevAltBase; }

  void SetRefHaplotype(std::string_view haplotype) { refHaplotype = Kmer::CanonicalSequence(haplotype); }
  void SetAltHaplotype(std::string_view haplotype) { altHaplotype = Kmer::CanonicalSequence(haplotype); }

  [[nodiscard]] auto ReferenceHaplotype() const -> std::string_view { return refHaplotype; }
  [[nodiscard]] auto AlternateHaplotype() const -> std::string_view { return altHaplotype; }

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

  // Canonical versions of the Reference and Alternate
  // haplotypes of atleast `k` length with Left & Right flank
  std::string refHaplotype;
  std::string altHaplotype;
};
}  // namespace lancet2
