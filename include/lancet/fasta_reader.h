#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "absl/status/statusor.h"
#include "lancet/contig_info.h"
#include "lancet/genomic_region.h"

namespace lancet {
class FastaReader {
 public:
  explicit FastaReader(const std::filesystem::path& ref);
  FastaReader() = delete;
  ~FastaReader();

  FastaReader(FastaReader&&) noexcept = default;
  auto operator=(FastaReader&&) noexcept -> FastaReader& = default;
  FastaReader(const FastaReader&) = delete;
  auto operator=(const FastaReader&) -> FastaReader& = delete;

  [[nodiscard]] auto RegionSequence(const GenomicRegion& region) const -> absl::StatusOr<std::string>;
  [[nodiscard]] auto ContigSequence(const std::string& contig) const -> absl::StatusOr<std::string>;

  [[nodiscard]] auto ContigsInfo() const -> std::vector<ContigInfo>;

  [[nodiscard]] auto ContigIDs() const -> absl::flat_hash_map<std::string, std::int64_t>;
  [[nodiscard]] auto ContigID(std::string_view contig) const -> absl::StatusOr<std::int64_t>;
  [[nodiscard]] auto ContigLength(std::string_view contig) const -> absl::StatusOr<std::int64_t>;

 private:
  class Impl;
  std::unique_ptr<Impl> pimpl;

  void Open(const std::filesystem::path& ref);
};
}  // namespace lancet
