#pragma once

#include <string>

#include "lancet2/sized_ints.h"

namespace lancet2 {
enum class Strand : bool { FWD, REV };
enum class SampleLabel : bool { NORMAL, TUMOR };
enum class KmerLabel : u8 { REFERENCE, NORMAL, TUMOR };
enum class EdgeKind : u8 { FF, FR, RF, RR };
enum class GraphEnd : bool { SOURCE, SINK };
enum class BuddyPosition : bool { FRONT, BACK };
enum class Allele : bool { REF, ALT };
enum class Haplotype : u8 { UNASSIGNED, FIRST, SECOND };
enum class TranscriptCode : u8 { REF_MATCH, SNV, INSERTION, DELETION, COMPLEX };
enum class VariantState : u8 { NONE, SOMATIC, NORMAL, SHARED };

[[nodiscard]] auto MakeEdgeKind(Strand first, Strand second) -> EdgeKind;
[[nodiscard]] auto SourceStrand(EdgeKind ek) -> Strand;
[[nodiscard]] auto DestStrand(EdgeKind ek) -> Strand;

[[nodiscard]] auto ReverseStrand(Strand s) -> Strand;
[[nodiscard]] auto ReverseSourceStrand(EdgeKind ek) -> EdgeKind;
[[nodiscard]] auto ReverseEdgeKind(EdgeKind ek) -> EdgeKind;

[[nodiscard]] auto ToString(Strand s) -> std::string;
[[nodiscard]] auto ToString(SampleLabel sl) -> std::string;
[[nodiscard]] auto ToString(KmerLabel kl) -> std::string;
[[nodiscard]] auto ToString(EdgeKind ek) -> std::string;
[[nodiscard]] auto ToString(GraphEnd ge) -> std::string;
[[nodiscard]] auto ToString(Allele al) -> std::string;
[[nodiscard]] auto ToString(Haplotype hp) -> std::string;
[[nodiscard]] auto ToString(TranscriptCode code) -> std::string;
[[nodiscard]] auto ToString(VariantState state) -> std::string;

extern const usize MOCK_SOURCE_ID;
extern const usize MOCK_SINK_ID;
}  // namespace lancet2
