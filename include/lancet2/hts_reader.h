#pragma once

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/types/span.h"
#include "lancet2/contig_info.h"
#include "lancet2/genomic_region.h"
#include "lancet2/hts_alignment.h"

namespace lancet2 {
class HtsReader {
 public:
  HtsReader(const std::filesystem::path& inpath, const std::filesystem::path& ref);
  ~HtsReader();
  HtsReader(HtsReader&&) noexcept = default;
  auto operator=(HtsReader&&) noexcept -> HtsReader& = default;

  HtsReader() = delete;
  HtsReader(const HtsReader&) = delete;
  auto operator=(const HtsReader&) -> HtsReader& = delete;

  auto JumpToContig(const std::string& contig) -> absl::Status;
  auto JumpToRegion(const GenomicRegion& region) -> absl::Status;
  auto SetBatchRegions(absl::Span<const GenomicRegion> regions) -> absl::Status;

  void ResetIterator();

  enum class IteratorState : int { INVALID = -2, DONE = -1, VALID = 0 };
  [[nodiscard]] auto GetNextAlignment(HtsAlignment* result, absl::Span<const std::string> fill_tags) -> IteratorState;

  [[nodiscard]] auto GetSampleNames() const -> std::vector<std::string>;
  [[nodiscard]] auto GetContigs() const -> std::vector<ContigInfo>;

  [[nodiscard]] auto GetContigIndex(const std::string& contig) const -> int;

 private:
  class Impl;
  std::unique_ptr<Impl> pimpl;
};

[[nodiscard]] auto TagPeekCheck(const std::filesystem::path& inpath, const std::filesystem::path& ref, const char* tag,
                                int max_alignments_to_read = 1000) -> bool;
}  // namespace lancet2
