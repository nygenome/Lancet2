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

  [[nodiscard]] auto chromosome() const noexcept -> std::string { return chromName; }
  [[nodiscard]] auto position() const noexcept -> std::size_t { return genomeRefPos; }
  [[nodiscard]] auto ref_start_offset() const noexcept -> std::size_t { return idxs.refStart; }
  [[nodiscard]] auto alt_start_offset() const noexcept -> std::size_t { return idxs.altStart; }

  [[nodiscard]] auto has_alt_coverage() const -> bool;

  [[nodiscard]] auto ref_end_offset() const noexcept -> std::size_t { return idxs.refEnd; }
  auto set_ref_end_offset(std::size_t val) -> Transcript&;

  [[nodiscard]] auto alt_end_offset() const noexcept -> std::size_t { return idxs.altEnd; }
  auto set_alt_end_offset(std::size_t val) -> Transcript&;

  [[nodiscard]] auto code() const noexcept -> TranscriptCode { return kind; }
  auto set_code(TranscriptCode val) -> Transcript&;

  [[nodiscard]] auto variant_coverage(SampleLabel label) const -> VariantHpCov;
  auto add_coverage(SampleLabel label, Allele al, const BaseHpCov& c) -> Transcript&;

  [[nodiscard]] auto str_result() const noexcept -> std::string;
  auto add_str_result(const TandemRepeatResult& val) -> Transcript&;

  [[nodiscard]] auto ref_seq() const noexcept -> std::string { return refSeq; }
  auto add_ref_base(const char& b) -> Transcript&;

  [[nodiscard]] auto alt_seq() const noexcept -> std::string { return altSeq; }
  auto add_alt_base(const char& b) -> Transcript&;

  [[nodiscard]] auto prev_ref_base() const noexcept -> char { return prevRefBase; }
  [[nodiscard]] auto prev_alt_base() const noexcept -> char { return prevAltBase; }

  [[nodiscard]] auto is_somatic() const noexcept -> bool { return isSomatic; }
  auto set_somatic_status(bool val = true) -> Transcript&;

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
