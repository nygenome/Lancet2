#pragma once

#include <string>
#include <string_view>

#include "lancet2/core_enums.h"
#include "lancet2/sized_ints.h"

namespace lancet2 {
class Kmer {
 public:
  Kmer(std::string_view sv) noexcept;  // NOLINT
  Kmer() = default;

  [[nodiscard]] auto CanMergeKmers(const Kmer& buddy, BuddyPosition merge_dir, bool reverse_buddy, usize k) const
      -> bool;

  void MergeBuddy(const Kmer& buddy, BuddyPosition dir, bool reverse_buddy, usize k);

  [[nodiscard]] auto GetFwdSeq() const -> std::string;
  [[nodiscard]] auto GetSeqView() const -> std::string_view { return seq; }
  [[nodiscard]] auto GetOrientation() const -> Strand { return strand; }
  [[nodiscard]] auto GetLength() const -> usize { return seq.length(); }
  [[nodiscard]] auto GetSize() const -> usize { return seq.size(); }
  [[nodiscard]] auto IsEmpty() const -> bool { return seq.empty(); }

  [[nodiscard]] auto GetHash() const -> u64;
  friend auto operator==(const Kmer& lhs, const Kmer& rhs) -> bool { return lhs.GetHash() == rhs.GetHash(); }
  friend auto operator!=(const Kmer& lhs, const Kmer& rhs) -> bool { return !(lhs == rhs); }

  [[nodiscard]] static auto CanonicalSequence(std::string_view sv) -> std::string;
  [[nodiscard]] static auto CanonicalSeqHash(std::string_view sv) -> u64;

 private:
  std::string seq;
  Strand strand = Strand::FWD;

  void Canonicalize(std::string_view sv);
};
}  // namespace lancet2
