#include "lancet2/kmer.h"

#include "absl/hash/internal/city.h"
#include "lancet2/merge_node_info.h"
#include "lancet2/utils.h"

namespace lancet2 {
Kmer::Kmer(std::string_view sv) noexcept { Canonicalize(sv); }

auto Kmer::CanMergeKmers(const Kmer& buddy, BuddyPosition merge_dir, bool reverse_buddy, usize k) const -> bool {
  return CanMergeSeqs(seq, buddy.GetSeqView(), merge_dir, reverse_buddy, k);
}

void Kmer::MergeBuddy(const Kmer& buddy, BuddyPosition dir, bool reverse_buddy, usize k) {
  seq.reserve(seq.length() + buddy.GetLength() - k + 1);
  MergeKmerSeqs(&seq, buddy.seq, dir, reverse_buddy, k);
}

auto Kmer::GetFwdSeq() const -> std::string { return strand == Strand::REV ? utils::RevComp(seq) : seq; }

void Kmer::Canonicalize(std::string_view sv) {
  const auto revComp = utils::RevComp(sv);
  if (sv < revComp) {
    seq = std::string(sv);
    strand = Strand::FWD;
    return;
  }

  seq = revComp;
  strand = Strand::REV;
}

auto Kmer::GetHash() const -> u64 {
  return absl::hash_internal::CityHash64WithSeeds(seq.c_str(), seq.length(), utils::PRIME0, utils::PRIME1);  // NOLINT
}

auto Kmer::CanonicalSequence(std::string_view sv) -> std::string {
  auto revComp = utils::RevComp(sv);
  if (sv < revComp) {
    return std::string(sv);
  }

  return revComp;
}

auto Kmer::CanonicalSeqHash(std::string_view sv) -> u64 {
  const auto cseq = CanonicalSequence(sv);
  return absl::hash_internal::CityHash64WithSeeds(cseq.c_str(), cseq.length(), utils::PRIME0, utils::PRIME1);  // NOLINT
}
}  // namespace lancet2
