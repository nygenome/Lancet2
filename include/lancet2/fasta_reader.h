#pragma once

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "absl/status/statusor.h"
#include "lancet2/contig_info.h"
#include "lancet2/genomic_region.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
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

  [[nodiscard]] auto ContigIDs() const -> absl::flat_hash_map<std::string, i64>;
  [[nodiscard]] auto ContigID(std::string_view contig) const -> absl::StatusOr<i64>;
  [[nodiscard]] auto ContigLength(std::string_view contig) const -> absl::StatusOr<i64>;

 private:
  class Impl;
  std::unique_ptr<Impl> pimpl;

  void Open(const std::filesystem::path& ref);
};
}  // namespace lancet2
