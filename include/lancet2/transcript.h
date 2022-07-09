#pragma once

#include <array>
#include <string>

#include "lancet2/core_enums.h"
#include "lancet2/sample_cov.h"
#include "lancet2/sized_ints.h"
#include "lancet2/tandem_repeat.h"
#include "lancet2/variant_hpcov.h"

namespace lancet2 {
/// 0-based offsets within reference and path sequences of the transcript
struct TranscriptOffsets {
  usize refStart = 0;
  usize altStart = 0;
  usize refEnd = 0;
  usize altEnd = 0;
};

struct TranscriptBases {
  char refBase = '!';
  char altBase = '!';
  char prevRefBase = '!';
  char prevAltBase = '!';
};

class Transcript {
 public:
  Transcript(std::string chrom, usize genome_ref_pos, TranscriptCode k, TranscriptOffsets offs, TranscriptBases bases,
             std::array<SampleCov, 2> covs, bool somatic_status);

  Transcript() = delete;

  [[nodiscard]] auto ChromName() const noexcept -> std::string { return chromName; }
  [[nodiscard]] auto Position() const noexcept -> usize { return genomeRefPos; }
  [[nodiscard]] auto RefStartOffset() const noexcept -> usize { return idxs.refStart; }
  [[nodiscard]] auto AltStartOffset() const noexcept -> usize { return idxs.altStart; }

  [[nodiscard]] auto HasAltCov() const -> bool;

  [[nodiscard]] auto RefEndOffset() const noexcept -> usize { return idxs.refEnd; }
  auto SetRefEndOffset(usize val) -> Transcript&;

  [[nodiscard]] auto AltEndOffset() const noexcept -> usize { return idxs.altEnd; }
  auto SetAltEndOffset(usize val) -> Transcript&;

  [[nodiscard]] auto Code() const noexcept -> TranscriptCode { return kind; }
  auto SetCode(TranscriptCode val) -> Transcript&;

  [[nodiscard]] auto VariantCov(SampleLabel label) const -> VariantHpCov;
  auto AddCov(SampleLabel label, Allele al, const BaseHpCov& c) -> Transcript&;

  [[nodiscard]] auto STRResult() const noexcept -> std::string;
  auto AddSTRResult(const TandemRepeatResult& val) -> Transcript&;

  [[nodiscard]] auto RefSeq() const noexcept -> std::string { return refSeq; }
  auto AddRefBase(const char& b) -> Transcript&;

  [[nodiscard]] auto AltSeq() const noexcept -> std::string { return altSeq; }
  auto AddAltBase(const char& b) -> Transcript&;

  [[nodiscard]] auto PrevRefBase() const noexcept -> char { return prevRefBase; }
  [[nodiscard]] auto PrevAltBase() const noexcept -> char { return prevAltBase; }

  [[nodiscard]] auto IsSomatic() const noexcept -> bool { return isSomatic; }
  auto SetSomaticStatus(bool val = true) -> Transcript&;

  [[nodiscard]] auto ComputeState() const -> VariantState;

 private:
  std::string chromName;
  usize genomeRefPos = 0;  // 1-based genome position for VCF
  TranscriptCode kind = TranscriptCode::REF_MATCH;
  TranscriptOffsets idxs;
  std::array<SampleCov, 2> sampleCovs;
  TandemRepeatResult strQry;
  std::string refSeq;
  std::string altSeq;
  char prevRefBase = 'N';
  char prevAltBase = 'N';
  bool isSomatic = false;
};
}  // namespace lancet2
