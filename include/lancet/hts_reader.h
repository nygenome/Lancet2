#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/types/span.h"
#include "lancet/contig_info.h"
#include "lancet/genomic_region.h"
#include "lancet/hts_alignment.h"
#include "lancet/statusor.h"

namespace lancet {
class HtsReader {
 public:
  HtsReader(const std::filesystem::path& inpath, const std::filesystem::path& ref);
  ~HtsReader();
  HtsReader(HtsReader&&) noexcept = default;
  auto operator=(HtsReader&&) noexcept -> HtsReader& = default;

  HtsReader() = delete;
  HtsReader(const HtsReader&) = delete;
  auto operator=(const HtsReader&) -> HtsReader& = delete;

  auto SetRegion(const std::string& contig) -> absl::Status;
  auto SetRegion(const GenomicRegion& region) -> absl::Status;
  auto SetRegions(absl::Span<const GenomicRegion> regions) -> absl::Status;

  void ResetIterator();

  enum class IteratorState : int { INVALID = -2, DONE = -1, VALID = 0 };
  [[nodiscard]] auto NextAlignment(HtsAlignment* result, absl::Span<const std::string> fill_tags) -> IteratorState;

  [[nodiscard]] auto SampleNames() const -> std::vector<std::string>;
  [[nodiscard]] auto ContigsInfo() const -> std::vector<ContigInfo>;

 private:
  class Impl;
  std::unique_ptr<Impl> pimpl;
};

[[nodiscard]] auto HasTag(const std::filesystem::path& inpath, const std::filesystem::path& ref, const char* tag,
                          int max_alignments_to_read = 1000) -> bool;
}  // namespace lancet
