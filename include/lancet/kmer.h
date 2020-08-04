#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

#include "lancet/core_enums.h"

namespace lancet {
class Kmer {
 public:
  Kmer(std::string_view sv) noexcept;  // NOLINT
  Kmer() = default;

  [[nodiscard]] auto CanMergeKmers(const Kmer& buddy, BuddyPosition merge_dir, bool reverse_buddy, std::size_t k) const
      -> bool;

  void MergeBuddy(const Kmer& buddy, BuddyPosition dir, bool reverse_buddy, std::size_t k);

  [[nodiscard]] auto FwdSeq() const -> std::string;
  [[nodiscard]] auto SeqView() const -> std::string_view { return seq; }
  [[nodiscard]] auto Orientation() const -> Strand { return strand; }
  [[nodiscard]] auto Length() const -> std::size_t { return seq.length(); }
  [[nodiscard]] auto Size() const -> std::size_t { return seq.size(); }
  [[nodiscard]] auto IsEmpty() const -> bool { return seq.empty(); }

  [[nodiscard]] auto ID() const -> std::uint64_t;
  friend auto operator==(const Kmer& lhs, const Kmer& rhs) -> bool { return lhs.ID() == rhs.ID(); }
  friend auto operator!=(const Kmer& lhs, const Kmer& rhs) -> bool { return !(lhs == rhs); }

  [[nodiscard]] static auto IsCanonical(std::string_view sv) -> bool;

 private:
  std::string seq;
  Strand strand = Strand::FWD;

  void Canonicalize(std::string_view sv);
};
}  // namespace lancet
