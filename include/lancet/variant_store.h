#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/synchronization/mutex.h"
#include "absl/types/span.h"
#include "lancet/cli_params.h"
#include "lancet/variant.h"
#include "lancet/vcf_writer.h"

namespace lancet {
class VariantStore {
 public:
  using WindowIds = std::vector<std::uint64_t>;
  using Container = absl::flat_hash_map<std::uint64_t, Variant>;

  explicit VariantStore(std::size_t num_windows, std::shared_ptr<const CliParams> p);
  VariantStore() = delete;

  [[nodiscard]] static auto BuildVcfHeader(const std::vector<std::string>& sample_names, const CliParams& p)
      -> std::string;

  auto AddVariant(std::size_t window_index, Variant&& variant) -> bool;

  auto FlushWindow(std::size_t window_index, VcfWriter* out,
                   const absl::flat_hash_map<std::string, std::int64_t>& contig_ids) -> bool;

  auto FlushAll(VcfWriter* out, const absl::flat_hash_map<std::string, std::int64_t>& contig_ids) -> bool;

 private:
  Container data{};
  absl::Mutex mutex;
  std::vector<WindowIds> windowIds;
  std::shared_ptr<const CliParams> params = nullptr;

  [[nodiscard]] auto FlushVariants(absl::Span<const std::uint64_t> variant_ids, VcfWriter* out,
                                   const absl::flat_hash_map<std::string, std::int64_t>& contig_ids) -> bool;
};
}  // namespace lancet
