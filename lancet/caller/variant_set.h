#ifndef SRC_LANCET_CALLER_VARIANT_SET_H_
#define SRC_LANCET_CALLER_VARIANT_SET_H_

#include <array>
#include <memory>
#include <string_view>
#include <vector>

#include "absl/container/btree_set.h"
#include "lancet/base/types.h"
#include "lancet/caller/msa_builder.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/core/window.h"

namespace lancet::caller {

class VariantSet {
 public:
  explicit VariantSet(const MsaBuilder& bldr, const core::Window& win);

  using Btree = absl::btree_set<RawVariant>;

  [[nodiscard]] auto begin() -> Btree::iterator { return mResultVariants.begin(); }
  [[nodiscard]] auto begin() const -> Btree::const_iterator { return mResultVariants.begin(); }
  [[nodiscard]] auto cbegin() const -> Btree::const_iterator { return mResultVariants.cbegin(); }

  [[nodiscard]] auto end() -> Btree::iterator { return mResultVariants.end(); }
  [[nodiscard]] auto end() const -> Btree::const_iterator { return mResultVariants.end(); }
  [[nodiscard]] auto cend() const -> Btree::const_iterator { return mResultVariants.cend(); }

  [[nodiscard]] auto rbegin() -> Btree::reverse_iterator { return mResultVariants.rbegin(); }
  [[nodiscard]] auto rbegin() const -> Btree::const_reverse_iterator { return mResultVariants.rbegin(); }
  [[nodiscard]] auto crbegin() const -> Btree::const_reverse_iterator { return mResultVariants.crbegin(); }

  [[nodiscard]] auto rend() -> Btree::reverse_iterator { return mResultVariants.rend(); }
  [[nodiscard]] auto rend() const -> Btree::const_reverse_iterator { return mResultVariants.rend(); }
  [[nodiscard]] auto crend() const -> Btree::const_reverse_iterator { return mResultVariants.crend(); }

  [[nodiscard]] auto IsEmpty() const -> bool { return mResultVariants.empty(); }
  [[nodiscard]] auto Count() const -> usize { return mResultVariants.size(); }

  using StartEndRange = std::array<i32, 2>;
  using VariantPtrs = std::vector<const RawVariant*>;
  [[nodiscard]] auto FindOverlappingVariants(usize hap_idx, const StartEndRange& interval) const -> VariantPtrs;

 private:
  absl::btree_set<RawVariant> mResultVariants;

  using StartEndIndices = std::array<usize, 2>;
  using VariationRanges = std::vector<StartEndIndices>;
  using PairwiseAlignment = std::array<std::string_view, 2>;
  [[nodiscard]] static auto FindVariationRanges(const PairwiseAlignment& alignment_view) -> VariationRanges;
  [[nodiscard]] static auto EndGapsFreeRange(const PairwiseAlignment& alignment_view) -> StartEndIndices;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_SET_H_
