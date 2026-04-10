#ifndef SRC_LANCET_CALLER_RAW_VARIANT_H_
#define SRC_LANCET_CALLER_RAW_VARIANT_H_

#include "lancet/base/sequence_complexity.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/graph_complexity.h"

#include "absl/container/flat_hash_map.h"
#include "absl/strings/str_cat.h"

#include <string>
#include <utility>
#include <vector>

#include <cstdint>

namespace lancet::caller {

class RawVariant {
 public:
  RawVariant() = default;

  enum class Type : i8 { REF = -1, SNV = 0, INS = 1, DEL = 2, MNP = 3, CPX = 4 };
  enum class State : i8 { NONE = -1, SHARED = 0, NORMAL = 1, TUMOR = 2, UNKNOWN = 3 };

  // ===========================================================================
  // STRICT SEQUENCE CORE CLASSIFICATION
  // ---------------------------------------------------------------------------
  // Extracts the mathematically pure Mutation Core by squeezing matching
  // 5' prefixes and 3' suffixes. Decouples Variant classification entirely
  // from arbitrary padding requirements and bounds logic heuristically.
  // ===========================================================================
  [[nodiscard]] static auto ClassifyVariant(std::string_view ref, std::string_view alt) -> Type;

  // ===========================================================================
  // DATA STRUCTURE CHOICE: `AltAllele` SUB-PAYLOAD
  // ---------------------------------------------------------------------------
  // Biologically, a single genomic locus can mutate into numerous varying alternative
  // forms across differing haplotypes or heterogeneous tumor populations (e.g., A -> C
  // in haplotype 1, but A -> T in haplotype 2).
  //
  // By creating an `AltAllele` sub-struct, we tightly pack every distinct mutated
  // alternative directly adjacent to each other internally within a parent variant block.
  // This maps 1:1 with the VCF spec behavior (where the ALT column holds "C,T"). It
  // prevents multiple artificially broken-apart biallelic records from misaligning or
  // receiving conflicting right-trim parsimony lengths later in the pipeline.
  // ===========================================================================
  struct AltAllele {
    std::string mSequence;

    // 3. MULTI-ALLELIC LOCAL MATRIX MAP (ALTs):
    // A single multi-allelic locus inherently possesses different string offsets
    // depending unilaterally on which specific structural path (Haplotype ID) a read traversed
    // to reach it. A 100bp insertion earlier in Haplotype 3 shifts this ALTs local index
    // mathematically by +100 relative to Haplotype 1!
    // -> Maps natively: Haplotype ID -> variant's exact Local Matrix Start on THAT string.
    absl::flat_hash_map<usize, usize>
        mLocalHapStart0Idxs;  // Moved up to satisfy 8B -> 4B -> 2B -> 1B constraints

    i64 mLength = -1;
    Type mType = Type::REF;

    // DATA STRUCTURE CHOICE: `absl::flat_hash_map`
    // We map a supporting physical haplotype ID -> to its 0-indexed start position within
    // the variant. Abseil's FlatHashMap utilizes Google's 'SwissTable' SIMD-accelerated
    // implementation. It stores keys/values contiguously in memory (Open Addressing), which
    // drastically outperforms std::map (O(logN) pointer chasing) and std::unordered_map
    // (cache-thrashing chunked linked lists). Max throughput here is essential.

    friend auto operator==(AltAllele const& lhs, AltAllele const& rhs) -> bool {
      return lhs.mSequence == rhs.mSequence && lhs.mLength == rhs.mLength && lhs.mType == rhs.mType;
    }

    friend auto operator<(AltAllele const& lhs, AltAllele const& rhs) -> bool {
      // Lexicographical sorting on `mSequence` establishes a strict deterministic
      // weak ordering. This guarantees deterministic binary behaviors when the parent
      // array of multiple ALTs is subsequently sorted in standard sets.
      return lhs.mSequence < rhs.mSequence;
    }

    template <typename HashState>
    friend auto AbslHashValue(HashState hash_state, AltAllele const& alt) -> HashState {
      // Leverage Abseil's strong native hashing composition frameworks
      return HashState::combine(std::move(hash_state), alt.mSequence, alt.mLength,
                                static_cast<i8>(alt.mType));
    }
  };

  usize mChromIndex = SIZE_MAX;

  // 1. GLOBAL GENOMIC COORDINATE: Strictly tracks where the variant structurally lands
  // on the actual reference genome. Only strictly used for VCF sorting and emitting.
  usize mGenomeChromPos1 = SIZE_MAX;

  // 2. LOCAL MATRIX COORDINATE (REF): The exact 0-indexed position within the
  // specific Reference string array spanning this exact Micro-Assembly window.
  // Necessary to exactly bind structural CIGAR strings backwards.
  usize mLocalRefStart0Idx = SIZE_MAX;

  std::string mChromName;
  // Contains the universal left-aligned bounding sequence encompassing ALL ALTs
  std::string mRefAllele;

  // DATA STRUCTURE CHOICE: `std::vector<AltAllele>`
  // Replaces the rigid single scalar ALTs. Vectors are cache optimized, avoiding heap
  // allocations for very small multiallelic blocks (as most sites have 1 or 2 alts max).
  std::vector<AltAllele> mAlts;

  // ── Graph complexity metrics (ML-ready orthogonal features) ────────────
  // Populated when --enable-graph-complexity-features is set.
  // 3 fields matching the GRAPH_CX VCF INFO tag.
  //
  // Coverage stability: GEI uses CovCV (σ/μ, self-normalizing ratio),
  // TipToPathCovRatio is a coverage ratio, MaxSingleDirDegree is pure
  // topology. All three are coverage-stable above 20×.
  //
  // Raw topology metrics (CC, BP, EdgeDensity, UnitigRatio, CoverageCv) are
  // mathematically collinear and compressed into the Graph Entanglement Index.
  // Color-based metrics (UnsharedColorRatio, ColorDiscordantBranches) are not
  // topological and are captured by other biologically relevant annotations.
  struct GraphMetrics {
    f64 mGraphEntanglementIndex = 0.0;  ///< GEI: log₁₀(1 + CC×BP×CovCV / UnitigRatio)
    f64 mTipToPathCovRatio =
        0.0;  ///< assembly tearing: tip cov / unitig cov (ratio, self-normalizing)
    usize mMaxSingleDirDegree =
        0;  ///< hub k-mer detection: max outgoing edges (topology, invariant)

    /// Format as 3 comma-separated values for VCF GRAPH_CX INFO tag.
    [[nodiscard]] auto FormatVcfValue() const -> std::string {
      return absl::StrCat(base::FormatComplexityScore(mGraphEntanglementIndex), ",",
                          base::FormatComplexityScore(mTipToPathCovRatio), ",",
                          mMaxSingleDirDegree);
    }
  };

  // ── Annotation fields (mutable — populated post-construction) ─────────
  // These do not participate in btree_set ordering and are annotated after
  // variant discovery, so they are mutable to allow modification through
  // const btree_set iterators without const_cast.

  mutable GraphMetrics mGraphMetrics;

  // ── Sequence complexity (11 ML-ready features) ────────────────────────
  // Populated when --enable-sequence-complexity-features is set.
  // Distilled from raw multi-scale metrics (HRun, Entropy, LongdustQ, TR motifs)
  // into 3 groups: Context (REF), Delta (ALT−REF), TR Motif (ALT).
  mutable base::SequenceComplexity mSeqCx;

  template <typename HashState>
  friend auto AbslHashValue(HashState hash_state, RawVariant const& var) -> HashState {
    return HashState::combine(std::move(hash_state), var.mChromIndex, var.mGenomeChromPos1,
                              var.mChromName, var.mRefAllele, var.mAlts);
  }

  friend auto operator==(RawVariant const& lhs, RawVariant const& rhs) -> bool {
    // Array equality operator `==` handles validating entire Multiallelic arrays instantly natively
    return lhs.mChromIndex == rhs.mChromIndex &&
           lhs.mGenomeChromPos1 == rhs.mGenomeChromPos1 &&
           lhs.mChromName == rhs.mChromName &&
           lhs.mRefAllele == rhs.mRefAllele &&
           lhs.mAlts == rhs.mAlts;
  }

  friend auto operator<(RawVariant const& lhs, RawVariant const& rhs) -> bool {
    if (lhs.mChromIndex != rhs.mChromIndex) return lhs.mChromIndex < rhs.mChromIndex;

    if (lhs.mGenomeChromPos1 != rhs.mGenomeChromPos1) {
      return lhs.mGenomeChromPos1 < rhs.mGenomeChromPos1;
    }

    if (lhs.mRefAllele != rhs.mRefAllele) return lhs.mRefAllele < rhs.mRefAllele;
    // Lexical vector comparison behaves safely
    return lhs.mAlts < rhs.mAlts;
  }

  friend auto operator!=(RawVariant const& lhs, RawVariant const& rhs) -> bool {
    return !(rhs == lhs);
  }
  friend auto operator>(RawVariant const& lhs, RawVariant const& rhs) -> bool { return rhs < lhs; }
  friend auto operator<=(RawVariant const& lhs, RawVariant const& rhs) -> bool {
    return !(rhs < lhs);
  }
  friend auto operator>=(RawVariant const& lhs, RawVariant const& rhs) -> bool {
    return !(lhs < rhs);
  }
};

// =========================================================================================
// STRICT SEQUENCE CORE CLASSIFICATION ENGINE
// =========================================================================================
// WHY DO WE SQUEEZE THE STRINGS AGAIN IF THE `VariantExtractor` FSM ALREADY APPLIES PARSIMONY?
//
// While `VariantBubble::NormalizeVcfParsimony` performs a simultaneous VT-style right-trim
// and left-alignment, it mathematically MUST execute identically universally across ALL alleles
// in a single Multi-Allelic Bubble block. Therefore, the global trimmer is inherently restricted
// by the structural topological bounds of the WIDEST allele in the cluster.
//
// -> THE MULTI-ALLELIC SHIELDING PROBLEM:
// Consider a bubble traversing 3 paths:
//    [REF]     :  A T G C
//    [ALT1]    :  A T
//    [ALT2]    :  A G G C
//
// Let's trace Global VCF Parsimony on ALT1 vs REF:
// It WANTS to right-trim `G` and `C` to isolate the `GC` deletion completely (ATGC -> AT).
// However, `ALT2` completely blocks the `C` and `G` from being universally erased, because
// `ALT2` actively utilizes them for its own structural integrity.
//
// Thus, VCF Parsimony halts prematurely for ALT1!
// When ALT1 evaluates its payload conventionally: `REF="ATGC"`, `ALT="AT"`.
// If we naively utilized length logic (`diff = -2` and `length > 1`), we would classify this
// erroneously as a `CPX` (Complex) mutation because the length boundary is artificially inflated!
//
// -> THE SEQUENCE CORE SOLUTION:
// By aggressively symmetrically squeezing matching 5' prefixes and 3' suffixes exclusively
// between the 1-on-1 pairs immediately prior to classification, we computationally decouple
// the biological Core from the VCF-Padding constraints.
//
//    [REF]     :  A T (G C)  ---> "GC"
//    [ALT1]    :  A T        ---> ""     ===> Result: Pure `DEL`!
//
// This executes in ultra-fast O(N) operations utilizing forward and reverse bounds pointers
// ensuring zero memory allocation overhead while achieving absolute biological mapping precision.
// =========================================================================================
inline auto RawVariant::ClassifyVariant(std::string_view ref, std::string_view alt) -> Type {
  usize start_match = 0;
  while (start_match < ref.length() &&
         start_match < alt.length() &&
         ref[start_match] == alt[start_match]) {
    start_match++;
  }

  if (start_match == ref.length() && start_match == alt.length()) {
    return Type::REF;
  }

  usize end_match =
      0;  // Ensures bounding pointers do not violently overlap or double-count characters
  while (end_match < (ref.length() - start_match) &&
         end_match < (alt.length() - start_match) &&
         ref[ref.length() - 1 - end_match] == alt[alt.length() - 1 - end_match]) {
    end_match++;
  }

  auto const ref_core_len = ref.length() - start_match - end_match;
  auto const alt_core_len = alt.length() - start_match - end_match;

  if (ref_core_len == 0 && alt_core_len > 0) {
    return Type::INS;
  }
  if (ref_core_len > 0 && alt_core_len == 0) {
    return Type::DEL;
  }
  if (ref_core_len > 0 && alt_core_len > 0) {
    if (ref_core_len == alt_core_len) {
      return (ref_core_len == 1) ? Type::SNV : Type::MNP;
    }
    return Type::CPX;  // Mixed/Complex concurrent insertion and deletion architectures
  }
  return Type::REF;
}

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_RAW_VARIANT_H_
