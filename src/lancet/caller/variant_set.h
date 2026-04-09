#ifndef SRC_LANCET_CALLER_VARIANT_SET_H_
#define SRC_LANCET_CALLER_VARIANT_SET_H_

#include <array>
#include <string_view>
#include <vector>

#include "absl/container/btree_map.h"
#include "absl/container/btree_set.h"
#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/caller/msa_builder.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/core/window.h"

namespace lancet::caller {

// ============================================================================
// VariantSet: Multiallelic Graph Extraction Engine
//
// Extracts biologically complex multiallelic variants natively from SPOA
// directed acyclic graphs by aggressively sweeping the topology concurrently.
// It tracks divergent structural paths per haplotype and mathematically merges 
// them into fully bundled `RawVariant` outputs with zero overlapping biases.
// ============================================================================
class VariantSet {
 public:
  VariantSet() = default;
  void ExtractVariantsFromGraph(const spoa::Graph& graph, const core::Window& win, usize ref_anchor_start);

  using BTree = absl::btree_set<RawVariant>;

  [[nodiscard]] auto begin() -> BTree::iterator { return mResultVariants.begin(); }
  [[nodiscard]] auto begin() const -> BTree::const_iterator { return mResultVariants.begin(); }
  [[nodiscard]] auto cbegin() const -> BTree::const_iterator { return mResultVariants.cbegin(); }

  [[nodiscard]] auto end() -> BTree::iterator { return mResultVariants.end(); }
  [[nodiscard]] auto end() const -> BTree::const_iterator { return mResultVariants.end(); }
  [[nodiscard]] auto cend() const -> BTree::const_iterator { return mResultVariants.cend(); }

  [[nodiscard]] auto rbegin() -> BTree::reverse_iterator { return mResultVariants.rbegin(); }
  [[nodiscard]] auto rbegin() const -> BTree::const_reverse_iterator { return mResultVariants.rbegin(); }
  [[nodiscard]] auto crbegin() const -> BTree::const_reverse_iterator { return mResultVariants.crbegin(); }

  [[nodiscard]] auto rend() -> BTree::reverse_iterator { return mResultVariants.rend(); }
  [[nodiscard]] auto rend() const -> BTree::const_reverse_iterator { return mResultVariants.rend(); }
  [[nodiscard]] auto crend() const -> BTree::const_reverse_iterator { return mResultVariants.crend(); }

  [[nodiscard]] auto IsEmpty() const -> bool { return mResultVariants.empty(); }
  [[nodiscard]] auto Count() const -> usize { return mResultVariants.size(); }

 private:
  absl::btree_set<RawVariant> mResultVariants;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_SET_H_
