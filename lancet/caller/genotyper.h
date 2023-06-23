#ifndef SRC_LANCET_CALLER_GENOTYPER_H_
#define SRC_LANCET_CALLER_GENOTYPER_H_

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/caller/variant_set.h"
#include "lancet/caller/variant_support.h"
#include "lancet/cbdg/read.h"
#include "minimap.h"

namespace lancet::caller {

class Genotyper {
 public:
  Genotyper();

  void SetNumSamples(const usize num_samples) { mNumSamples = num_samples; }

  using Reads = absl::Span<const cbdg::Read>;
  using RefHap = const hts::Reference::Region&;
  using AltHaps = absl::Span<const std::string>;

  using PerSampleVariantEvidence = absl::flat_hash_map<std::string_view, std::unique_ptr<VariantSupport>>;
  using Result = absl::flat_hash_map<const RawVariant*, PerSampleVariantEvidence>;
  [[nodiscard]] auto Genotype(RefHap item, AltHaps alts, Reads reads, const VariantSet& vset) -> Result;

  class AlnInfo {
   public:
    // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
    i32 mRefStart = -1;
    i32 mQryStart = -1;
    i32 mRefEnd = -1;
    i32 mQryEnd = -1;
    i32 mDpScore = -1;
    f64 mGcIden = 0.0;
    usize mHapIdx = 0;
    std::string mCsTag;
    // NOLINTEND(misc-non-private-member-variables-in-classes)

    [[nodiscard]] auto IsEmpty() const noexcept -> bool;

    using QryStartAllele = std::pair<usize, Allele>;
    using SupportsInfo = absl::flat_hash_map<const RawVariant*, QryStartAllele>;
    void AddSupportingInfo(SupportsInfo& supports, const VariantSet& vset) const;

   private:
    using StartEndIndices = std::array<usize, 2>;
    using IdentityRanges = std::vector<StartEndIndices>;
    using RefQryIdentityRanges = std::array<IdentityRanges, 2>;

    [[nodiscard]] auto FindIdentityRanges() const -> RefQryIdentityRanges;
    [[nodiscard]] static auto FindQueryStart(const RefQryIdentityRanges& ref_qry_equal_ranges,
                                             const StartEndIndices& allele_span) -> std::optional<usize>;
  };

 private:
  struct MmIdxDeleter {
    void operator()(mm_idx_t* idx) noexcept { mm_idx_destroy(idx); }
  };

  struct MmTbufDeleter {
    void operator()(mm_tbuf_t* tbuf) noexcept { mm_tbuf_destroy(tbuf); }
  };

  using MappingOpts = std::unique_ptr<mm_mapopt_t>;
  using IndexingOpts = std::unique_ptr<mm_idxopt_t>;
  using ThreadBuffer = std::unique_ptr<mm_tbuf_t, MmTbufDeleter>;
  using Minimap2Index = std::unique_ptr<mm_idx_t, MmIdxDeleter>;
  static constexpr usize REF_HAP_IDX = 0;

  usize mNumSamples = 0;
  std::vector<Minimap2Index> mIndices;
  MappingOpts mMappingOpts = std::make_unique<mm_mapopt_t>();
  IndexingOpts mIndexingOpts = std::make_unique<mm_idxopt_t>();
  ThreadBuffer mThreadBuffer = ThreadBuffer(mm_tbuf_init());

  void ResetData(RefHap ref, AltHaps alts);

  [[nodiscard]] auto AlignRead(const cbdg::Read& read) -> std::vector<AlnInfo>;

  using SupportsInfo = AlnInfo::SupportsInfo;
  static void AddToTable(Result& result, const cbdg::Read& read, const SupportsInfo& read_supports);
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_GENOTYPER_H_
