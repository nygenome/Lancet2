#ifndef SRC_LANCET_HTS_ALIGNMENT_H_
#define SRC_LANCET_HTS_ALIGNMENT_H_

#include <string>
#include <string_view>
#include <vector>

extern "C" {
#include "htslib/hts.h"
#include "htslib/sam.h"
}

#include "absl/container/flat_hash_set.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/hts/cigar_unit.h"
#include "lancet/hts/reference.h"
#include "spdlog/fmt/bundled/core.h"

namespace lancet::hts {

/// Zero-copy, lightweight proxy over a `bam1_t*` record managed by the HTS iterator.
///
/// IMPORTANT LIFETIME WARNING:
/// This object holds a non-owning pointer (`mRawAln`) to the `bam1_t` block that is structurally
/// managed by the `hts::Iterator`. The pointer becomes INVALID after any subsequent `++itr` call
/// on the parent iterator. Consequently:
///   - Do NOT store `Alignment` objects beyond the current loop iteration.
///   - Do NOT hold references/pointers to data returned by `QnameView()`, `CigarData()`, etc.
///     across iterator increments.
///   - Any data that must outlive the iterator step (sequence, qualities) must be explicitly
///     extracted via `BuildSequence()` / `BuildQualities()` which perform deep copies on demand.
///
/// In practice, all existing call sites consume `Alignment` within range-for loop bodies and
/// immediately decompose relevant fields into owned types (e.g. `cbdg::Read`), so this zero-copy
/// approach is safe for the current codebase.
class Alignment {
 public:
  enum class Fields : u16 {
    CORE_QNAME = SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_RNEXT | SAM_PNEXT | SAM_TLEN,
    SEQ_QUAL = CORE_QNAME | SAM_SEQ | SAM_QUAL,
    CIGAR_SEQ_QUAL = SEQ_QUAL | SAM_CIGAR,
    AUX_RGAUX = CIGAR_SEQ_QUAL | SAM_AUX | SAM_RGAUX,
  };

  enum class Strand : bool { FWD = true, REV = false };

  class BitwiseFlag {
   public:
    BitwiseFlag(u16 flags) : mFlag(flags) {}
    BitwiseFlag() = default;

    [[nodiscard]] auto GetStrand() const noexcept -> Strand;
    [[nodiscard]] auto GetMateStrand() const noexcept -> Strand;
    [[nodiscard]] auto IsFwdStrand() const noexcept -> bool;
    [[nodiscard]] auto IsRevStrand() const noexcept -> bool;
    [[nodiscard]] auto IsMateFwdStrand() const noexcept -> bool;
    [[nodiscard]] auto IsMateRevStrand() const noexcept -> bool;
    [[nodiscard]] auto IsQcFail() const noexcept -> bool;
    [[nodiscard]] auto IsDuplicate() const noexcept -> bool;
    [[nodiscard]] auto IsPrimary() const noexcept -> bool;
    [[nodiscard]] auto IsSecondary() const noexcept -> bool;
    [[nodiscard]] auto IsSupplementary() const noexcept -> bool;
    [[nodiscard]] auto IsMapped() const noexcept -> bool;
    [[nodiscard]] auto IsUnmapped() const noexcept -> bool;
    [[nodiscard]] auto IsMateMapped() const noexcept -> bool;
    [[nodiscard]] auto IsMateUnmapped() const noexcept -> bool;
    [[nodiscard]] auto IsPairedInSequencing() const noexcept -> bool;
    [[nodiscard]] auto IsMappedProperPair() const noexcept -> bool;
    [[nodiscard]] auto IsRead1() const noexcept -> bool;
    [[nodiscard]] auto IsRead2() const noexcept -> bool;
    [[nodiscard]] auto HasFlagsSet(u16 check_flags) const noexcept -> bool;
    [[nodiscard]] auto HasFlagsUnset(u16 check_flags) const noexcept -> bool;

    [[nodiscard]] operator u16() const noexcept { return mFlag; }

    auto operator==(const BitwiseFlag& rhs) const -> bool { return mFlag == rhs.mFlag; }
    auto operator!=(const BitwiseFlag& rhs) const -> bool { return !(rhs == *this); }

   private:
    u16 mFlag = 0;
  };

  // ---------------------------------------------------------------------------
  // Core fields: these read directly from cached scalar fields (always populated)
  // ---------------------------------------------------------------------------
  [[nodiscard]] auto StartPos0() const noexcept -> i64 { return mStart0; }
  [[nodiscard]] auto MateStartPos0() const noexcept -> i64 { return mMateStart0; }
  [[nodiscard]] auto InsertSize() const noexcept -> i64 { return mInsertSize; }
  [[nodiscard]] auto ChromIndex() const noexcept -> i32 { return mChromIdx; }
  [[nodiscard]] auto MateChromIndex() const noexcept -> i32 { return mMateChromIdx; }
  [[nodiscard]] auto Flag() const noexcept -> BitwiseFlag { return mSamFlag; }
  [[nodiscard]] auto FlagRaw() const noexcept -> u16 { return mSamFlag; }
  [[nodiscard]] auto MapQual() const noexcept -> u8 { return mMapQual; }

  // ---------------------------------------------------------------------------
  // Zero-copy proxies: route directly through the underlying bam1_t* payload.
  // The returned views are only valid while this Alignment (and its backing
  // iterator) remain at the current position.
  // ---------------------------------------------------------------------------

  /// Zero-copy view of the query name directly from the bam1_t record.
  [[nodiscard]] auto QnameView() const noexcept -> std::string_view;

  /// Returns the query sequence length from the bam1_t core fields.
  [[nodiscard]] auto Length() const noexcept -> usize;

  /// Zero-copy view of the raw CIGAR array from the bam1_t record.
  [[nodiscard]] auto CigarData() const -> std::vector<CigarUnit>;

  [[nodiscard]] auto CigarString() const -> std::string;

  // ---------------------------------------------------------------------------
  // On-demand deep extraction: these perform allocations and must be used when
  // the data needs to outlive the current iterator position.
  // ---------------------------------------------------------------------------

  /// Decodes the 4-bit packed BAM sequence into a full ASCII string.
  /// This allocates a new std::string on every call.
  [[nodiscard]] auto BuildSequence() const -> std::string;

  /// Copies the raw quality values into a new vector.
  /// This allocates a new std::vector<u8> on every call.
  [[nodiscard]] auto BuildQualities() const -> std::vector<u8>;

  struct MateInfo {
    i32 mChromIndex = -1;
    i64 mMateStartPos0 = -1;
  };

  [[nodiscard]] auto MateLocation() const noexcept -> MateInfo;
  [[nodiscard]] auto MateOverlapsRegion(const Reference::Region& region) const noexcept -> bool;

  [[nodiscard]] auto OverlapsRegion(const Reference::Region& region) const noexcept -> bool;

  [[nodiscard]] auto IsEmpty() const noexcept -> bool;

  // ---------------------------------------------------------------------------
  // Aux tag access: routes directly through bam_aux_get() on the raw record,
  // bypassing the old cached vector of AuxTag objects.
  // ---------------------------------------------------------------------------

  /// Check if a two-character auxiliary tag exists in this alignment's raw data.
  [[nodiscard]] auto HasTag(std::string_view tag_name) const noexcept -> bool;

  /// Retrieve the value of an auxiliary tag. Supported types:
  ///   - i64 for integer tags (c/C/s/S/i/I)
  ///   - f64 for float tags (f/d)
  ///   - std::string_view for string tags (Z/H) -- valid only while iterator is at current position
  template <typename TagResultValue>
  [[nodiscard]] auto GetTag(std::string_view tag_name) const -> absl::StatusOr<TagResultValue> {
    const auto* raw_aux = FindRawTag(tag_name);
    if (raw_aux == nullptr) {
      const auto msg = fmt::format("Tag {} is not present in the alignment record", tag_name);
      return absl::Status(absl::StatusCode::kNotFound, msg);
    }
    return ExtractTagValue<TagResultValue>(raw_aux, tag_name);
  }

  [[nodiscard]] auto GetSoftClips(std::vector<u32>* clip_sizes, std::vector<u32>* read_positions,
                                  std::vector<u32>* genome_positions, bool use_padded = false) const -> bool;

  [[nodiscard]] auto ToString(const Reference& ref) const -> std::string;

  auto operator==(const Alignment& rhs) const -> bool;
  auto operator!=(const Alignment& rhs) const -> bool;

 private:
  // Cached scalar fields from bam1_t::core (always populated, cheap to copy)
  i64 mStart0 = -1;
  i64 mMateStart0 = -1;
  i64 mInsertSize = -1;
  i32 mChromIdx = -1;
  i32 mMateChromIdx = -1;
  u16 mSamFlag = 0;
  u8 mMapQual = 0;

  /// Non-owning pointer to the bam1_t block managed by the Iterator/Extractor.
  /// WARNING: This pointer is invalidated on the next iterator increment (++itr).
  /// See the class-level documentation for full lifetime semantics.
  bam1_t* mRawAln = nullptr;

  friend class Iterator;
  using TagNamesSet = absl::flat_hash_set<std::string>;

  Alignment() = default;

  void ClearAllFields();
  void PopulateFromRaw(bam1_t* aln);

  /// Direct bam_aux_get lookup. Returns nullptr if tag is not present.
  [[nodiscard]] auto FindRawTag(std::string_view tag_name) const noexcept -> const u8*;

  /// Type-dispatch helper for extracting a typed value from a raw aux pointer.
  template <typename TagResultValue>
  [[nodiscard]] static auto ExtractTagValue(const u8* raw_aux, std::string_view tag_name)
      -> absl::StatusOr<TagResultValue>;
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_ALIGNMENT_H_
