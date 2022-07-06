#pragma once

#include <string>
#include <string_view>
#include <utility>

#include "absl/container/flat_hash_map.h"
#include "absl/status/statusor.h"
#include "lancet2/cigar.h"
#include "lancet2/core_enums.h"
#include "lancet2/genomic_region.h"
#include "lancet2/read_info.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class HtsAlignment {
 public:
  HtsAlignment() = default;

  [[nodiscard]] auto ReadName() const -> std::string { return readName; }
  [[nodiscard]] auto ContigName() const -> std::string { return contig; }
  [[nodiscard]] auto MateContigName() const -> std::string { return mateContig; }
  [[nodiscard]] auto ReadSequence() const -> std::string { return readSequence; }
  [[nodiscard]] auto ReadQuality() const -> std::string { return readQuality; }

  [[nodiscard]] auto StartPosition0() const -> i64 { return startPosition0; }
  [[nodiscard]] auto EndPosition0() const -> i64 { return endPosition0; }
  [[nodiscard]] auto MateStartPosition0() const -> i64 { return mateStartPosition0; }
  [[nodiscard]] auto MappingQuality() const -> u8 { return mappingQuality; }

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

  [[nodiscard]] auto Length() const -> usize { return readSequence.length(); }

  // NOTE: Returned tag data only valid until HtsReader's alignment data is valid
  [[nodiscard]] auto HasTag(std::string_view tag) const -> bool { return tagsData.contains(tag); }
  [[nodiscard]] auto TagData(std::string_view tag) const -> absl::StatusOr<const u8*> {
    const auto itr = tagsData.find(tag);
    if (itr == tagsData.end()) return absl::NotFoundError("requested tag is not present");
    return itr->second;
  }

  [[nodiscard]] auto BuildReadInfo(SampleLabel label, u8 min_bq, u8 max_kmer_size) -> ReadInfo;

  [[nodiscard]] auto SoftClips(std::vector<u32>* clip_sizes, std::vector<u32>* read_positions,
                               std::vector<u32>* genome_positions, bool use_padded = false) -> bool;

  void SetReadName(std::string rname) { readName = std::move(rname); }
  void SetContig(std::string ctg) { contig = std::move(ctg); }
  void SetMateContig(std::string matectg) { mateContig = std::move(matectg); }
  void SetReadSequence(std::string seq) { readSequence = std::move(seq); }
  void SetReadQuality(std::string qual) { readQuality = std::move(qual); }

  void SetStartPosition0(i64 start) { startPosition0 = start; }
  void SetEndPosition0(i64 end) { endPosition0 = end; }
  void SetMateStartPosition0(i64 matestart) { mateStartPosition0 = matestart; }
  void SetMappingQuality(u8 mapqual) { mappingQuality = mapqual; }

  void SetSamFlags(u16 flags) { samFlags = flags; }
  void SetCigar(AlignmentCigar cig) { cigar = std::move(cig); }
  void SetTagData(std::string_view tag, const u8* data) { tagsData[std::string(tag)] = data; }

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

  i64 startPosition0 = -1;
  i64 endPosition0 = -1;
  i64 mateStartPosition0 = -1;
  u8 mappingQuality = 0;
  u16 samFlags = 0;

  AlignmentCigar cigar;
  absl::flat_hash_map<std::string, const u8*> tagsData;
};
}  // namespace lancet2
