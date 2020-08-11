#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include "absl/base/thread_annotations.h"
#include "absl/container/flat_hash_map.h"
#include "absl/synchronization/mutex.h"
#include "absl/types/span.h"
#include "lancet/cli_params.h"
#include "lancet/variant.h"

namespace lancet {
class VariantStore {
 public:
  using WindowIds = std::vector<std::uint64_t>;

  explicit VariantStore(std::size_t num_windows, std::shared_ptr<const CliParams> p);
  VariantStore() = delete;

  [[nodiscard]] static auto BuildVcfHeader(const std::vector<std::string>& sample_names, const CliParams& p)
      -> std::string;

  auto ABSL_LOCKS_EXCLUDED(mutex) AddVariant(std::size_t window_index, Variant&& variant) -> bool;

  auto ABSL_LOCKS_EXCLUDED(mutex) FlushWindow(std::size_t window_index, std::ostream& out,
                                              const absl::flat_hash_map<std::string, std::int64_t>& ctg_ids) -> bool;

  auto ABSL_LOCKS_EXCLUDED(mutex)
      FlushAll(std::ostream& out, const absl::flat_hash_map<std::string, std::int64_t>& ctg_ids) -> bool;

 private:
  mutable absl::Mutex mutex;
  std::vector<WindowIds> windowIds ABSL_GUARDED_BY(mutex);
  absl::flat_hash_map<std::uint64_t, Variant> data ABSL_GUARDED_BY(mutex);
  std::shared_ptr<const CliParams> params = nullptr;

  [[nodiscard]] ABSL_EXCLUSIVE_LOCKS_REQUIRED(mutex) auto FlushVariants(
      absl::Span<const std::uint64_t> variant_ids, std::ostream& out,
      const absl::flat_hash_map<std::string, std::int64_t>& ctg_ids) -> bool;

  [[nodiscard]] static auto IsVariantLesser(const Variant& v1, const Variant& v2,
                                            const absl::flat_hash_map<std::string, std::int64_t>& ctg_ids) -> bool;
};
}  // namespace lancet
