#include "lancet/cbdg/kmer.h"

#include <array>

#include "absl/strings/string_view.h"
#include "lancet/base/hash.h"
#include "lancet/base/rev_comp.h"

namespace {

// Get overlapping prefix and suffix portions of adjacent kmers
inline auto OvlPrefix(std::string_view data, const usize kval) { return data.substr(0, kval - 1); }
inline auto OvlSuffix(std::string_view data, const usize kval) { return data.substr(data.size() - kval + 1, kval - 1); }

// Get non-overlapping prefix and suffix portions of adjacent kmers
inline auto NonOvlPrefix(std::string_view data, const usize kval) { return data.substr(0, data.size() - kval + 1); }
inline auto NonOvlSuffix(std::string_view data, const usize kval) {
  return data.substr(kval - 1, data.size() - kval + 1);
}

/// Logic to merge sequence and quality data from k2 into k1. See comments below for details
/// https://github.com/GATB/bcalm/blob/v2.2.3/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void MergeCords(std::string& k1_dflt, std::string_view k2_dflt, const lancet::cbdg::EdgeKind ekind, const usize currk) {
  k1_dflt.reserve(k1_dflt.length() + k2_dflt.length() - currk);
  using lancet::cbdg::EdgeKind;
  switch (ekind) {
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
      k1_dflt.append(NonOvlSuffix(k2_dflt, currk));
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
      k1_dflt.append(NonOvlSuffix(RevComp(k2_dflt), currk));
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
      k1_dflt.insert(0, NonOvlPrefix(RevComp(k2_dflt), currk));
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
      k1_dflt.insert(0, NonOvlPrefix(k2_dflt, currk));
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
      mIdentifier = HashStr64(seq);
      mDfltSeq = seq;
      break;

    case Sign::MINUS:
      mIdentifier = HashStr64(rc_seq);
      mDfltSeq = std::move(rc_seq);
      break;
  }
}

void Kmer::Merge(const Kmer& other, const EdgeKind conn_kind, usize currk) {
  if (IsEmpty()) {
    mDfltSign = other.mDfltSign;
    mIdentifier = other.mIdentifier;
    mDfltSeq = other.mDfltSeq;
    return;
  }

  MergeCords(mDfltSeq, other.mDfltSeq, conn_kind, currk);
}

auto Kmer::SignFor(const Ordering order) const noexcept -> Sign {
  return order == Ordering::DEFAULT ? mDfltSign : RevSign(mDfltSign);
}

auto Kmer::SequenceFor(const Ordering order) const -> std::string {
  return order == Ordering::DEFAULT ? mDfltSeq : RevComp(mDfltSeq);
}

}  // namespace lancet::cbdg
