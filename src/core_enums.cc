#include "lancet2/core_enums.h"

#include "absl/hash/internal/city.h"
#include "lancet2/utils.h"

namespace lancet2 {
auto MakeEdgeKind(Strand first, Strand second) -> EdgeKind {
  if (first == Strand::FWD && second == Strand::FWD) return EdgeKind::FF;
  if (first == Strand::FWD && second == Strand::REV) return EdgeKind::FR;
  if (first == Strand::REV && second == Strand::FWD) return EdgeKind::RF;
  return EdgeKind::RR;
}

auto SourceStrand(EdgeKind ek) -> Strand {
  return (ek == EdgeKind::FF || ek == EdgeKind::FR) ? Strand::FWD : Strand::REV;
}

auto DestStrand(EdgeKind ek) -> Strand {
  return (ek == EdgeKind::FF || ek == EdgeKind::RF) ? Strand::FWD : Strand::REV;
}

auto ReverseStrand(Strand s) -> Strand { return s == Strand::FWD ? Strand::REV : Strand::FWD; }

auto ReverseSourceStrand(EdgeKind ek) -> EdgeKind {
  switch (ek) {
    case EdgeKind::FF:
      return EdgeKind::RF;
    case EdgeKind::FR:
      return EdgeKind::RR;
    case EdgeKind::RF:
      return EdgeKind::FF;
    case EdgeKind::RR:
    default:
      return EdgeKind::FR;
  }
}

auto ReverseEdgeKind(EdgeKind ek) -> EdgeKind {
  if (ek == EdgeKind::FR || ek == EdgeKind::RF) return ek;
  return ek == EdgeKind::FF ? EdgeKind::RR : EdgeKind::FF;
}

auto ToString(Strand s) -> std::string { return s == Strand::FWD ? "F" : "R"; }
auto ToString(SampleLabel sl) -> std::string { return sl == SampleLabel::NORMAL ? "NML" : "TMR"; }

auto ToString(KmerLabel kl) -> std::string {
  switch (kl) {
    case KmerLabel::REFERENCE:
      return "ref";
    case KmerLabel::NORMAL:
      return "nml";
    case KmerLabel::TUMOR:
    default:
      return "tmr";
  }
}

auto ToString(EdgeKind ek) -> std::string {
  switch (ek) {
    case EdgeKind::FF:
      return "FF";
    case EdgeKind::FR:
      return "FR";
    case EdgeKind::RF:
      return "RF";
    case EdgeKind::RR:
    default:
      return "RR";
  }
}

auto ToString(GraphEnd ge) -> std::string { return ge == GraphEnd::SOURCE ? "source" : "sink"; }
auto ToString(Allele al) -> std::string { return al == Allele::REF ? "REF" : "ALT"; }

auto ToString(Haplotype hp) -> std::string {
  switch (hp) {
    case Haplotype::FIRST:
      return "1";
    case Haplotype::SECOND:
      return "2";
    case Haplotype::UNASSIGNED:
    default:
      return "0";
  }
}

auto ToString(TranscriptCode code) -> std::string {
  switch (code) {
    case TranscriptCode::SNV:
      return "snv";
    case TranscriptCode::INSERTION:
      return "ins";
    case TranscriptCode::DELETION:
      return "del";
    case TranscriptCode::COMPLEX:
      return "complex";
    case TranscriptCode::REF_MATCH:
    default:
      return "ref_match";
  }
}

auto ToString(VariantState state) -> std::string {
  switch (state) {
    case VariantState::SOMATIC:
      return "SOMATIC";
    case VariantState::NORMAL:
      return "NORMAL";
    case VariantState::SHARED:
      return "SHARED";
    case VariantState::NONE:
    default:
      return "NONE";
  }
}

const std::size_t MOCK_SOURCE_ID =
    absl::hash_internal::CityHash64WithSeeds("MOCK_SOURCE", 11, utils::PRIME_0, utils::PRIME_1);  // NOLINT

const std::size_t MOCK_SINK_ID =
    absl::hash_internal::CityHash64WithSeeds("MOCK_SINK", 9, utils::PRIME_0, utils::PRIME_1);  // NOLINT
}  // namespace lancet2
