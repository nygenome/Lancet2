#include "lancet/cbdg/kmer.h"

#include <array>

#include "lancet/base/assert.h"
#include "lancet/base/hash.h"
#include "lancet/base/rev_comp.h"

namespace {

// Get overlapping prefix and suffix portions of adjacent kmers
inline auto OvlPrefix(std::string_view data, const usize kval) { return data.substr(0, kval - 1); }
inline auto OvlSuffix(std::string_view data, const usize kval) { return data.substr(data.size() - kval + 1, kval - 1); }

// Get non-overlapping prefix and suffix portions of adjacent kmers
inline auto NonOvlPrefix(const absl::Cord& data, const usize kval) { return data.Subcord(0, data.size() - kval + 1); }
inline auto NonOvlSuffix(const absl::Cord& data, const usize kval) {
  return data.Subcord(kval - 1, data.size() - kval + 1);
}

/// Logic to merge sequence and quality data from k2 into k1. See comments below for details
/// https://github.com/GATB/bcalm/blob/v2.2.3/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
void MergeCords(absl::Cord& k1_dflt, absl::Cord& k1_oppo, const absl::Cord& k2_dflt, const absl::Cord& k2_oppo,
                const lancet::cbdg::EdgeKind edge_kind, const usize currk) {
#ifdef LANCET_DEVELOP_MODE
  const auto str_k1_dflt = std::string(k1_dflt);
  const auto str_k1_oppo = std::string(k1_oppo);
  const auto str_k2_dflt = std::string(k2_dflt);
  const auto str_k2_oppo = std::string(k2_oppo);

  LANCET_ASSERT(!str_k1_dflt.empty())
  LANCET_ASSERT(!str_k2_dflt.empty())
  LANCET_ASSERT(!str_k1_oppo.empty())
  LANCET_ASSERT(!str_k2_oppo.empty())

  std::string_view first_lhs;
  std::string_view first_rhs;
  std::string_view second_lhs;
  std::string_view second_rhs;
#endif

  using lancet::cbdg::EdgeKind;
  switch (edge_kind) {
    case EdgeKind::PLUS_PLUS:
      // src_fwd  5' ACCCACCTAATCCGACGCCGGTGCACCCGGGATACCGCATCTGTCTACC 3'
      // src_rev  3' GGTAGACAGATGCGGTATCCCGGGTGCACCGGCGTCGGATTAGGTGGGT 5'
      //---------------------------------------------------------------------
      // k1-->k2
      // k1_dflt  5' ACCCACCTAATCCGACGCCGGTGCACCCGGGAT                 3'
      // k2_dflt  5'                  CCGGTGCACCCGGGATACCGCATCTGTCTACC 3'
      //---------------------------------------------------------------------
      // k2-->k1
      // k1_oppo  3'                 ATCCCGGGTGCACCGGCGTCGGATTAGGTGGGT 5'
      // k2_oppo  3' GGTAGACAGATGCGGTATCCCGGGTGCACCGG                  5'
#ifdef LANCET_DEVELOP_MODE
      first_lhs = OvlSuffix(str_k1_dflt, currk);
      first_rhs = OvlPrefix(str_k2_dflt, currk);
      second_lhs = OvlPrefix(str_k1_oppo, currk);
      second_rhs = OvlSuffix(str_k2_oppo, currk);

      LANCET_ASSERT(first_lhs == first_rhs)
      LANCET_ASSERT(second_lhs == second_rhs)
#endif

      k1_dflt.Append(NonOvlSuffix(k2_dflt, currk));
      k1_oppo.Prepend(NonOvlPrefix(k2_oppo, currk));
      break;

    case EdgeKind::PLUS_MINUS:
      // src_fwd  5' CCTTACGGGAATAGGTGTGCCCCAATTTCTCCCATGAGGGTAACCTCGT 3'
      // src_rev  3' ACGAGGTTACCCTCATGGGAGAAATTGGGGCACACCTATTCCCGTAAGG 5'
      //---------------------------------------------------------------------
      // k1-->k2
      // k1_dflt  5' CCTTACGGGAATAGGTGTGCCCCAATTTCTCCC                 3'
      // k2_oppo  5'            TAGGTGTGCCCCAATTTCTCCCATGAGGGTAACCTCGT 3'
      //---------------------------------------------------------------------
      // k2-->k1
      // k1_oppo  3'                 GGGAGAAATTGGGGCACACCTATTCCCGTAAGG 5'
      // k2_dflt  3' ACGAGGTTACCCTCATGGGAGAAATTGGGGCACACCTA            5'
#ifdef LANCET_DEVELOP_MODE
      first_lhs = OvlSuffix(str_k1_dflt, currk);
      first_rhs = OvlPrefix(str_k2_oppo, currk);
      second_lhs = OvlPrefix(str_k1_oppo, currk);
      second_rhs = OvlSuffix(str_k2_dflt, currk);

      LANCET_ASSERT(first_lhs == first_rhs)
      LANCET_ASSERT(second_lhs == second_rhs)
#endif

      k1_dflt.Append(NonOvlSuffix(k2_oppo, currk));
      k1_oppo.Prepend(NonOvlPrefix(k2_dflt, currk));
      break;

    case EdgeKind::MINUS_PLUS:
      // src_fwd  5' GGAACTTTTGTACTATGAATTACGTAAAAAAGGGCTTGTATAGAAATCG 3'
      // src_rev  3' CGATTTCTATACAAGCCCTTTTTTACGTAATTCATAGTACAAAAGTTCC 5'
      //---------------------------------------------------------------------
      // k1-->k2
      // k1_oppo  5' GGAACTTTTGTACTATGAATTACGTAAAAAAGG                 3'
      // k2_dflt  5'                  AATTACGTAAAAAAGGGCTTGTATAGAAATCG 3'
      //---------------------------------------------------------------------
      // k2-->k1
      // k1_dflt  3'                 CCTTTTTTACGTAATTCATAGTACAAAAGTTCC 5'
      // k2_oppo  3' CGATTTCTATACAAGCCCTTTTTTACGTAATT                  5'
#ifdef LANCET_DEVELOP_MODE
      first_lhs = OvlSuffix(str_k1_oppo, currk);
      first_rhs = OvlPrefix(str_k2_dflt, currk);
      second_lhs = OvlPrefix(str_k1_dflt, currk);
      second_rhs = OvlSuffix(str_k2_oppo, currk);

      LANCET_ASSERT(first_lhs == first_rhs)
      LANCET_ASSERT(second_lhs == second_rhs)
#endif

      k1_dflt.Prepend(NonOvlPrefix(k2_oppo, currk));
      k1_oppo.Append(NonOvlSuffix(k2_dflt, currk));
      break;

    case EdgeKind::MINUS_MINUS:
      // src_fwd  5' TCCGTATGTTGCAGCATTGTGCTACTCGTTTGGATACAGGTAAGGGCGT 3'
      // src_rev  3' ACGCCCTTACCTGTATCCAAACGAGTAGCACAATGCTGCAACATACGGA 5'
      //---------------------------------------------------------------------
      // k1<--k2
      // k1_oppo  5' TCCGTATGTTGCAGCATTGTGCTACTCGTTTGG                 3'
      // k2_oppo  5'                  TGTGCTACTCGTTTGGATACAGGTAAGGGCGT 3'
      //---------------------------------------------------------------------
      // k2<--k1
      // k1_dflt  3'                 CCAAACGAGTAGCACAATGCTGCAACATACGGA 5'
      // k2_dflt  3' ACGCCCTTACCTGTATCCAAACGAGTAGCACA                  5'
#ifdef LANCET_DEVELOP_MODE
      first_lhs = OvlSuffix(str_k1_oppo, currk);
      first_rhs = OvlPrefix(str_k2_oppo, currk);
      second_lhs = OvlPrefix(str_k1_dflt, currk);
      second_rhs = OvlSuffix(str_k2_dflt, currk);

      LANCET_ASSERT(first_lhs == first_rhs)
      LANCET_ASSERT(second_lhs == second_rhs)
#endif

      k1_dflt.Prepend(NonOvlPrefix(k2_dflt, currk));
      k1_oppo.Append(NonOvlSuffix(k2_oppo, currk));
      break;
  }
}

}  // namespace

namespace lancet::cbdg {

Kmer::Kmer(std::string_view seq) {
  auto rc_seq = RevComp(seq);
  mDfltSign = seq < rc_seq ? Sign::PLUS : Sign::MINUS;

  switch (mDfltSign) {
    case Sign::PLUS:
      mIdentifier = HashStr(seq);
      mDfltSeq = seq;
      mOppoSeq = std::move(rc_seq);
      break;

    case Sign::MINUS:
      mIdentifier = HashStr(rc_seq);
      mDfltSeq = std::move(rc_seq);
      mOppoSeq = seq;
      break;
  }
}

void Kmer::Merge(const Kmer& other, const EdgeKind conn_kind, usize currk) {
  if (IsEmpty()) {
    mDfltSign = other.mDfltSign;
    mIdentifier = other.mIdentifier;
    mDfltSeq = other.mDfltSeq;
    mOppoSeq = other.mOppoSeq;
    return;
  }

  MergeCords(mDfltSeq, mOppoSeq, other.mDfltSeq, other.mOppoSeq, conn_kind, currk);
}

auto Kmer::Length() const -> usize {
  LANCET_ASSERT(mDfltSeq.size() == mOppoSeq.size())
  return mDfltSeq.size();
}

auto Kmer::SignFor(const Ordering order) const noexcept -> Sign {
  return order == Ordering::DEFAULT ? mDfltSign : RevSign(mDfltSign);
}

auto Kmer::SequenceFor(const Ordering order) const -> std::string {
  return order == Ordering::DEFAULT ? static_cast<std::string>(mDfltSeq) : static_cast<std::string>(mOppoSeq);
}

auto Kmer::CordDataFor(const Ordering order) const -> absl::Cord {
  return order == Ordering::DEFAULT ? mDfltSeq : mOppoSeq;
}

[[nodiscard]] auto SlidingKmers(std::string_view seq, const usize window) -> absl::FixedArray<Kmer> {
  if (seq.length() < window) {
    return absl::FixedArray<Kmer>(0);
  }

  const auto end_position = seq.length() - window;
  absl::FixedArray<Kmer> result(end_position + 1);

  for (usize offset = 0; offset <= end_position; ++offset) {
    result[offset] = Kmer(absl::ClippedSubstr(seq, offset, window));
    LANCET_ASSERT(result[offset].Length() == window)
  }

  return result;
}

}  // namespace lancet::cbdg
