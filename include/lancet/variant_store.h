#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "absl/base/thread_annotations.h"
#include "absl/container/flat_hash_map.h"
#include "absl/synchronization/mutex.h"
#include "absl/types/span.h"
#include "lancet/cli_params.h"
#include "lancet/ref_window.h"
#include "lancet/variant.h"
#include "lancet/vcf_writer.h"

namespace lancet {
class VariantStore {
 public:
  // ChromosomeName -> ChromosomeIndex in reference FASTA
  using ContigIDs = absl::flat_hash_map<std::string, std::int64_t>;

  explicit VariantStore(std::shared_ptr<const CliParams> p);
  VariantStore() = delete;

  [[nodiscard]] static auto GetHeader(const std::vector<std::string>& sample_names, const CliParams& p) -> std::string;

  auto ABSL_LOCKS_EXCLUDED(mutex) AddVariant(Variant&& variant) -> bool;
  auto ABSL_LOCKS_EXCLUDED(mutex) FlushWindow(const RefWindow& w, VcfWriter* out, const ContigIDs& ctg_ids) -> bool;
  auto ABSL_LOCKS_EXCLUDED(mutex) FlushAll(VcfWriter* out, const ContigIDs& ctg_ids) -> bool;

 private:
  mutable absl::Mutex mutex;
  absl::flat_hash_map<VariantID, Variant> data ABSL_GUARDED_BY(mutex);
  std::shared_ptr<const CliParams> params = nullptr;

  [[nodiscard]] ABSL_EXCLUSIVE_LOCKS_REQUIRED(mutex) auto Flush(absl::Span<const VariantID> ids, VcfWriter* out,
                                                                const ContigIDs& ctg_ids) -> bool;

  [[nodiscard]] static auto IsVariant1LessThan2(const Variant& v1, const Variant& v2, const ContigIDs& ctg_ids) -> bool;
  [[nodiscard]] static auto IsVariantInOrBefore(const Variant& v, const RefWindow& w, const ContigIDs& ctg_ids) -> bool;
};
}  // namespace lancet
