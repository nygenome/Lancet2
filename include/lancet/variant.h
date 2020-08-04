#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>

#include "lancet/cli_params.h"
#include "lancet/transcript.h"
#include "lancet/variant_hpcov.h"

namespace lancet {
class Variant {
 public:
  explicit Variant(const Transcript& transcript, std::size_t kmer_size);
  Variant() = delete;

  [[nodiscard]] auto MakeVcfLine(const CliParams& params) const -> std::string;
  [[nodiscard]] auto ID() const -> std::uint64_t;
  [[nodiscard]] auto ComputeState() const -> VariantState;

  std::string ChromName;   // NOLINT
  std::size_t Position;    // NOLINT
  std::string RefAllele;   // NOLINT
  std::string AltAllele;   // NOLINT
  TranscriptCode Kind;     // NOLINT
  std::string STRResult;   // NOLINT
  std::size_t Length;      // NOLINT
  std::size_t KmerSize;    // NOLINT
  VariantHpCov TumorCov;   // NOLINT
  VariantHpCov NormalCov;  // NOLINT

  auto operator==(const Variant& other) const -> bool { return this->ID() == other.ID(); }
  auto operator!=(const Variant& other) const -> bool { return this->ID() != other.ID(); }
  auto operator<(const Variant& other) const -> bool;

 private:
  static auto Genotype(int ref, int alt) -> std::string;
  static auto BuildSampleFormat(const VariantHpCov& v, bool is_tenx_mode) -> std::string;
};
}  // namespace lancet
