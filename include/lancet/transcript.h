#pragma once

#include <array>
#include <cstddef>
#include <string>

#include "lancet/core_enums.h"
#include "lancet/sample_cov.h"
#include "lancet/tandem_repeat.h"
#include "lancet/variant_hpcov.h"

namespace lancet {
/// 0-based offsets within reference and path sequences of the transcript
struct TranscriptOffsets {
  std::size_t refStart = 0;
  std::size_t altStart = 0;
  std::size_t refEnd = 0;
  std::size_t altEnd = 0;
};

struct TranscriptBases {
  char refBase = '!';
  char altBase = '!';
  char prevRefBase = '!';
  char prevAltBase = '!';
};

class Transcript {
 public:
  Transcript(std::string chrom, std::size_t genome_ref_pos, TranscriptCode k, TranscriptOffsets offs,
             TranscriptBases bases, std::array<SampleCov, 2> covs, bool somatic_status);

  Transcript() = delete;

  [[nodiscard]] auto ChromName() const noexcept -> std::string { return chromName; }
  [[nodiscard]] auto Position() const noexcept -> std::size_t { return genomeRefPos; }
  [[nodiscard]] auto RefStartOffset() const noexcept -> std::size_t { return idxs.refStart; }
  [[nodiscard]] auto AltStartOffset() const noexcept -> std::size_t { return idxs.altStart; }

  [[nodiscard]] auto HasAltCov() const -> bool;

  [[nodiscard]] auto RefEndOffset() const noexcept -> std::size_t { return idxs.refEnd; }
  auto SetRefEndOffset(std::size_t val) -> Transcript&;

  [[nodiscard]] auto AltEndOffset() const noexcept -> std::size_t { return idxs.altEnd; }
  auto SetAltEndOffset(std::size_t val) -> Transcript&;

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

 private:
  std::string chromName;
  std::size_t genomeRefPos = 0;  // 1-based genome position for VCF
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
}  // namespace lancet
