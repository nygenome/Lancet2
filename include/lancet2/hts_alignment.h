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

  [[nodiscard]] auto GetReadName() const -> std::string { return readName; }
  [[nodiscard]] auto GetContigName() const -> std::string { return contig; }
  [[nodiscard]] auto GetMateContigName() const -> std::string { return mateContig; }
  [[nodiscard]] auto GetReadSequence() const -> std::string { return readSequence; }
  [[nodiscard]] auto GetReadQuality() const -> std::string { return readQuality; }

  [[nodiscard]] auto GetStartPos0() const -> i64 { return startPosition0; }
  [[nodiscard]] auto GetEndPos0() const -> i64 { return endPosition0; }
  [[nodiscard]] auto GetMateStartPos0() const -> i64 { return mateStartPosition0; }
  [[nodiscard]] auto GetMappingQual() const -> u8 { return mappingQuality; }

  [[nodiscard]] auto GetCigarData() const -> AlignmentCigar { return cigar; }
  [[nodiscard]] auto GetReadStrand() const -> Strand;

  [[nodiscard]] auto GetMateRegion() const -> GenomicRegion {
    return {mateContig, static_cast<u32>(mateStartPosition0 + 1), static_cast<u32>(mateStartPosition0 + 1)};
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
    return (startPosition0 >= (region.GetStartPos1() - 1)) && (endPosition0 <= (region.GetEndPos1() - 1));
  }

  [[nodiscard]] auto PercentOverlapWith(const GenomicRegion& region) const -> double;

  [[nodiscard]] auto GetLength() const -> usize { return readSequence.length(); }

  // NOTE: Returned tag data only valid until HtsReader's alignment data is valid
  [[nodiscard]] auto HasTag(std::string_view tag) const -> bool { return tagsData.contains(tag); }
  [[nodiscard]] auto GetTagData(std::string_view tag) const -> absl::StatusOr<const u8*> {
    const auto itr = tagsData.find(tag);
    if (itr == tagsData.end()) return absl::NotFoundError("requested tag is not present");
    return itr->second;
  }

  [[nodiscard]] auto BuildReadInfo(SampleLabel label, u8 min_bq, u8 max_kmer_size) const -> ReadInfo;

  [[nodiscard]] auto GetSoftClips(std::vector<u32>* clip_sizes, std::vector<u32>* read_positions,
                                  std::vector<u32>* genome_positions, bool use_padded = false) const -> bool;

  void SetReadName(std::string rname) { readName = std::move(rname); }
  void SetContigName(std::string ctg) { contig = std::move(ctg); }
  void SetMateContigName(std::string matectg) { mateContig = std::move(matectg); }
  void SetReadSequence(std::string seq) { readSequence = std::move(seq); }
  void SetReadQuality(std::string qual) { readQuality = std::move(qual); }

  void SetStartPos0(i64 start) { startPosition0 = start; }
  void SetEndPos0(i64 end) { endPosition0 = end; }
  void SetMateStartPos0(i64 matestart) { mateStartPosition0 = matestart; }
  void SetMappingQual(u8 mapqual) { mappingQuality = mapqual; }

  void SetSamFlags(u16 flags) { samFlags = flags; }
  void SetCigarData(AlignmentCigar cig) { cigar = std::move(cig); }
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
