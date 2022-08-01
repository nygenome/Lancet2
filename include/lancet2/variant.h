#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>

#include "lancet2/cli_params.h"
#include "lancet2/core_enums.h"
#include "lancet2/transcript.h"
#include "lancet2/variant_hpcov.h"

namespace lancet2 {
using VariantID = u64;

class Variant {
 public:
  Variant(const Transcript& transcript, usize kmer_size, VariantHpCov tmrCov, VariantHpCov nmlCov);
  Variant() = delete;

  [[nodiscard]] auto MakeVcfLine(const CliParams& params) const -> std::string;
  [[nodiscard]] auto ID() const -> VariantID;
  [[nodiscard]] auto ComputeState() const -> VariantState;

  std::string ChromName;   // NOLINT
  usize Position;          // NOLINT
  std::string RefAllele;   // NOLINT
  std::string AltAllele;   // NOLINT
  TranscriptCode Kind;     // NOLINT
  std::string STRResult;   // NOLINT
  usize Length;            // NOLINT
  usize KmerSize;          // NOLINT
  VariantHpCov TumorCov;   // NOLINT
  VariantHpCov NormalCov;  // NOLINT

  float TmrRefQual = 0.0F;  // NOLINT
  float TmrAltQual = 0.0F;  // NOLINT
  float NmlRefQual = 0.0F;  // NOLINT
  float NmlAltQual = 0.0F;  // NOLINT

  auto operator==(const Variant& other) const -> bool { return this->ID() == other.ID(); }
  auto operator!=(const Variant& other) const -> bool { return this->ID() != other.ID(); }

 private:
  static auto Genotype(int ref, int alt) -> std::string;
  static auto BuildSampleFormat(const VariantHpCov& v, bool is_tenx_mode) -> std::string;
};
}  // namespace lancet2
