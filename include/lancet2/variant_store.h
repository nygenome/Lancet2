#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include "absl/base/thread_annotations.h"
#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"
#include "lancet2/cli_params.h"
#include "lancet2/ref_window.h"
#include "lancet2/spinlock.h"
#include "lancet2/variant.h"

namespace lancet2 {
class VariantStore {
 public:
  // ChromosomeName -> ChromosomeIndex in reference FASTA
  using ContigIDs = absl::flat_hash_map<std::string, std::int64_t>;

  explicit VariantStore(std::shared_ptr<const CliParams> p);
  VariantStore() = delete;

  [[nodiscard]] static auto GetHeader(const std::vector<std::string>& sample_names, const CliParams& p) -> std::string;

  /// Tries to add variants to store if no other thread is writing to it and
  /// returns `true` iff variants were added and `false` otherwise
  [[nodiscard]] auto TryAddVariants(absl::Span<const Variant> variants) -> bool;

  /// Blocks if any other thread is writing to store and waits to write variants to store
  void ForceAddVariants(absl::Span<const Variant> variants);

  auto FlushWindow(const RefWindow& w, std::ostream& out, const ContigIDs& ctg_ids) -> bool;
  auto FlushAll(std::ostream& out, const ContigIDs& ctg_ids) -> bool;

 private:
  // https://nathanpetersen.com/2019/02/17/optimizing-for-workloads-linux-spinlocks-vs-mutexes/
  mutable utils::SpinLock spinLock = utils::SpinLock();
  absl::flat_hash_map<VariantID, Variant> data;
  std::shared_ptr<const CliParams> params = nullptr;

  [[nodiscard]] auto Flush(absl::Span<const VariantID> ids, std::ostream& out, const ContigIDs& ctg_ids) -> bool;

  [[nodiscard]] static auto IsVariant1LessThan2(const Variant& v1, const Variant& v2, const ContigIDs& ctg_ids) -> bool;
  [[nodiscard]] static auto IsVariantInOrBefore(const Variant& v, const RefWindow& w, const ContigIDs& ctg_ids) -> bool;

  void UnsafeAddVariantBatch(absl::Span<const Variant> variants);
};
}  // namespace lancet2
