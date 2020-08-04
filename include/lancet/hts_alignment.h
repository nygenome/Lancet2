#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <utility>

#include "absl/container/flat_hash_map.h"
#include "lancet/cigar.h"
#include "lancet/core_enums.h"
#include "lancet/genomic_region.h"
#include "lancet/read_info.h"
#include "lancet/statusor.h"

namespace lancet {
class HtsAlignment {
 public:
  HtsAlignment() = default;

  [[nodiscard]] auto ReadName() const -> std::string { return readName; }
  [[nodiscard]] auto ContigName() const -> std::string { return contig; }
  [[nodiscard]] auto MateContigName() const -> std::string { return mateContig; }
  [[nodiscard]] auto ReadSequence() const -> std::string { return readSequence; }
  [[nodiscard]] auto ReadQuality() const -> std::string { return readQuality; }

  [[nodiscard]] auto StartPosition0() const -> std::int64_t { return startPosition0; }
  [[nodiscard]] auto EndPosition0() const -> std::int64_t { return endPosition0; }
  [[nodiscard]] auto MateStartPosition0() const -> std::int64_t { return mateStartPosition0; }
  [[nodiscard]] auto MappingQuality() const -> std::uint8_t { return mappingQuality; }

  [[nodiscard]] auto CigarData() const -> AlignmentCigar { return cigar; }
  [[nodiscard]] auto ReadStrand() const -> Strand;

  [[nodiscard]] auto MateRegion() const -> GenomicRegion {
    return GenomicRegion(mateContig, mateStartPosition0 + 1, mateStartPosition0 + 1);
  }

  [[nodiscard]] auto IsDuplicate() const -> bool;
  [[nodiscard]] auto IsSupplementary() const -> bool;
  [[nodiscard]] auto IsPrimary() const -> bool;
  [[nodiscard]] auto IsSecondary() const -> bool;
  [[nodiscard]] auto IsQcFailed() const -> bool;
  [[nodiscard]] auto IsUnmapped() const -> bool;
  [[nodiscard]] auto IsMateUnmapped() const -> bool;
  [[nodiscard]] auto IsReverseStrand() const -> bool;
  [[nodiscard]] auto IsMateReverseStrand() const -> bool;
  [[nodiscard]] auto IsPaired() const -> bool;
  [[nodiscard]] auto IsProperPair() const -> bool;
  [[nodiscard]] auto IsRead1() const -> bool;
  [[nodiscard]] auto IsRead2() const -> bool;

  [[nodiscard]] auto IsWithinRegion(const GenomicRegion& region) const -> bool {
    // GenomicRegion is 1-based, Alignment is 0-based
    return (startPosition0 >= (region.StartPosition1() - 1)) && (endPosition0 <= (region.EndPosition1() - 1));
  }

  [[nodiscard]] auto Length() const -> std::size_t { return readSequence.length(); }

  // NOTE: Returned tag data only valid until HtsReader's alignment data is valid
  [[nodiscard]] auto HasTag(std::string_view tag) const -> bool { return tagsData.contains(tag); }
  [[nodiscard]] auto TagData(std::string_view tag) const -> StatusOr<const std::uint8_t*> {
    const auto itr = tagsData.find(tag);
    if (itr == tagsData.end()) return absl::NotFoundError("requested tag is not present");
    return itr->second;
  }

  [[nodiscard]] auto BuildReadInfo(SampleLabel label, std::uint8_t min_bq, std::uint8_t max_kmer_size) -> ReadInfo;

  [[nodiscard]] auto SoftClips(std::vector<std::uint32_t>* clip_sizes, std::vector<std::uint32_t>* read_positions,
                               std::vector<std::uint32_t>* genome_positions, bool use_padded = false) -> bool;

  void SetReadName(std::string rname) { readName = std::move(rname); }
  void SetContig(std::string ctg) { contig = std::move(ctg); }
  void SetMateContig(std::string matectg) { mateContig = std::move(matectg); }
  void SetReadSequence(std::string seq) { readSequence = std::move(seq); }
  void SetReadQuality(std::string qual) { readQuality = std::move(qual); }

  void SetStartPosition0(std::int64_t start) { startPosition0 = start; }
  void SetEndPosition0(std::int64_t end) { endPosition0 = end; }
  void SetMateStartPosition0(std::int64_t matestart) { mateStartPosition0 = matestart; }
  void SetMappingQuality(std::uint8_t mapqual) { mappingQuality = mapqual; }

  void SetSamFlags(std::uint16_t flags) { samFlags = flags; }
  void SetCigar(AlignmentCigar cig) { cigar = std::move(cig); }
  void SetTagData(std::string_view tag, const std::uint8_t* data) { tagsData[std::string(tag)] = data; }

  void Clear() {
    readName.erase(readName.begin(), readName.end());
    contig.erase(contig.begin(), contig.end());
    mateContig.erase(mateContig.begin(), mateContig.end());
    readSequence.erase(readSequence.begin(), readSequence.end());
    readQuality.erase(readQuality.begin(), readQuality.end());
    startPosition0 = -1;
    endPosition0 = -1;
    mateStartPosition0 = -1;
    mappingQuality = 0;
    samFlags = 0;
    cigar.clear();
    tagsData.erase(tagsData.begin(), tagsData.end());
  }

 private:
  std::string readName;
  std::string contig;
  std::string mateContig;
  std::string readSequence;
  std::string readQuality;

  std::int64_t startPosition0 = -1;
  std::int64_t endPosition0 = -1;
  std::int64_t mateStartPosition0 = -1;
  std::uint8_t mappingQuality = 0;
  std::uint16_t samFlags = 0;

  AlignmentCigar cigar;
  absl::flat_hash_map<std::string, const std::uint8_t*> tagsData;
};
}  // namespace lancet
