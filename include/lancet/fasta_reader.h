#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "lancet/contig_info.h"
#include "lancet/genomic_region.h"
#include "lancet/statusor.h"

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

  [[nodiscard]] auto RegionSequence(const GenomicRegion& region) const -> StatusOr<std::string>;
  [[nodiscard]] auto ContigSequence(const std::string& contig) const -> StatusOr<std::string>;

  [[nodiscard]] auto ContigsInfo() const -> std::vector<ContigInfo>;

  [[nodiscard]] auto ContigId(std::string_view contig) const -> StatusOr<std::int64_t>;
  [[nodiscard]] auto ContigLength(std::string_view contig) const -> StatusOr<std::int64_t>;

 private:
  class Impl;
  std::unique_ptr<Impl> pimpl;

  void Open(const std::filesystem::path& ref);
};
}  // namespace lancet
