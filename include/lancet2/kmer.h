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

  [[nodiscard]] auto FwdSeq() const -> std::string;
  [[nodiscard]] auto SeqView() const -> std::string_view { return seq; }
  [[nodiscard]] auto Orientation() const -> Strand { return strand; }
  [[nodiscard]] auto Length() const -> usize { return seq.length(); }
  [[nodiscard]] auto Size() const -> usize { return seq.size(); }
  [[nodiscard]] auto IsEmpty() const -> bool { return seq.empty(); }

  [[nodiscard]] auto ID() const -> u64;
  friend auto operator==(const Kmer& lhs, const Kmer& rhs) -> bool { return lhs.ID() == rhs.ID(); }
  friend auto operator!=(const Kmer& lhs, const Kmer& rhs) -> bool { return !(lhs == rhs); }

  [[nodiscard]] static auto IsCanonical(std::string_view sv) -> bool;

 private:
  std::string seq;
  Strand strand = Strand::FWD;

  void Canonicalize(std::string_view sv);
};
}  // namespace lancet2
