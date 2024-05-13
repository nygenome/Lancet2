#ifndef SRC_LANCET_CBDG_KMER_H_
#define SRC_LANCET_CBDG_KMER_H_

#include <array>
#include <string>
#include <string_view>

#include "lancet/base/types.h"

namespace lancet::cbdg {

enum class EdgeKind : u8 { PLUS_PLUS = 0, PLUS_MINUS = 1, MINUS_PLUS = 2, MINUS_MINUS = 3 };

class Kmer {
 public:
  // Default sequence is lexicographically lower than opposite sequence
  enum class Ordering : bool { DEFAULT = true, OPPOSITE = false };
  [[nodiscard]] static auto RevOrdering(const Ordering ord) -> Ordering {
    return ord == Ordering::DEFAULT ? Ordering::OPPOSITE : Ordering::DEFAULT;
  }

  //  PLUS -> default seq in original orientation of the source read
  // MINUS -> default seq in rev_comp orientation of the source read
  enum class Sign : bool { PLUS = true, MINUS = false };
  [[nodiscard]] static auto RevSign(const Sign sign) -> Sign { return sign == Sign::PLUS ? Sign::MINUS : Sign::PLUS; }

  Kmer() = default;
  explicit Kmer(std::string_view seq);

  void Merge(const Kmer& other, EdgeKind conn_kind, usize currk);

  [[nodiscard]] auto SignFor(Ordering order) const noexcept -> Sign;
  [[nodiscard]] auto SequenceFor(Ordering order) const -> std::string;

  [[nodiscard]] auto Identifier() const noexcept -> u64 { return mIdentifier; }
  [[nodiscard]] auto Length() const -> usize { return mDfltSeq.length(); }
  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return mDfltSeq.empty() && mIdentifier == 0; }
  friend auto operator==(const Kmer& lhs, const Kmer& rhs) noexcept -> bool { return lhs.mDfltSeq == rhs.mDfltSeq; }
  friend auto operator!=(const Kmer& lhs, const Kmer& rhs) noexcept -> bool { return !(rhs == lhs); }

 private:
  Sign mDfltSign = Sign::PLUS;
  u64 mIdentifier = 0;
  std::string mDfltSeq;
};

[[nodiscard]] static inline auto MakeFwdEdgeKind(const std::array<Kmer::Sign, 2>& sign_pair) -> EdgeKind {
  const auto [src_sign, dst_sign] = sign_pair;
  switch (src_sign) {
    case (Kmer::Sign::PLUS):
      return dst_sign == Kmer::Sign::PLUS ? EdgeKind::PLUS_PLUS : EdgeKind::PLUS_MINUS;
    case (Kmer::Sign::MINUS):
      return dst_sign == Kmer::Sign::PLUS ? EdgeKind::MINUS_PLUS : EdgeKind::MINUS_MINUS;
  }
}

[[nodiscard]] static inline auto SplitIntoSignPair(const EdgeKind kind) -> std::array<Kmer::Sign, 2> {
  switch (kind) {
    case EdgeKind::PLUS_PLUS:
      return {Kmer::Sign::PLUS, Kmer::Sign::PLUS};

    case EdgeKind::PLUS_MINUS:
      return {Kmer::Sign::PLUS, Kmer::Sign::MINUS};

    case EdgeKind::MINUS_PLUS:
      return {Kmer::Sign::MINUS, Kmer::Sign::PLUS};

    default:
      return {Kmer::Sign::MINUS, Kmer::Sign::MINUS};
  }
}

[[nodiscard]] static inline auto RevEdgeKind(EdgeKind kind) -> EdgeKind {
  switch (kind) {
    case EdgeKind::PLUS_PLUS:
      return EdgeKind::MINUS_MINUS;

    case EdgeKind::MINUS_MINUS:
      return EdgeKind::PLUS_PLUS;

    default:
      return kind;
  }
}

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_KMER_H_
