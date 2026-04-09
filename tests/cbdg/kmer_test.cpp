#include "lancet/cbdg/kmer.h"

#include "lancet/base/assert.h"
#include "lancet/base/rev_comp.h"
#include "lancet/base/sliding.h"
#include "lancet/base/types.h"

#include "absl/container/fixed_array.h"
#include "catch_amalgamated.hpp"

#include <absl/strings/string_view.h>
#include <array>
#include <random>
#include <ranges>
#include <string>
#include <string_view>

using namespace lancet::cbdg;

namespace {

inline auto GenerateRandomDnaSequence(usize const seq_len) -> std::string {
  static constexpr std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};

  std::random_device device;
  std::mt19937_64 generator(device());

  std::uniform_int_distribution<usize> base_chooser(0, 3);
  std::string result(seq_len, 'N');

  for (usize idx = 0; idx < seq_len; ++idx) {
    result[idx] = BASES.at(base_chooser(generator));
  }

  return result;
}

inline auto MatchesOneOfTwo(std::string_view result, std::array<std::string_view, 2> const& values)
    -> bool {
  return (result == values[0]) || (result == values[1]);
}

[[nodiscard]] inline auto SlidingKmers(std::string_view seq, usize const window)
    -> absl::FixedArray<Kmer> {
  if (seq.length() < window) {
    return absl::FixedArray<Kmer>(0);
  }

  auto const end_position = seq.length() - window;
  absl::FixedArray<Kmer> result(end_position + 1);

  for (usize offset = 0; offset <= end_position; ++offset) {
    result[offset] = Kmer(absl::ClippedSubstr(seq, offset, window));
    LANCET_ASSERT(result[offset].Length() == window)
  }

  return result;
}

}  // namespace

static constexpr auto NUM_RANDOM_ITERATIONS = 100;
static constexpr auto DFLT_ORD = Kmer::Ordering::DEFAULT;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Can merge two adjacent equal sized kmers", "[lancet][cbdg][Kmer]") {
  static constexpr auto SEQ_LEN = 151;
  static constexpr auto KMER_SIZE = 11;

  absl::FixedArray<std::string, NUM_RANDOM_ITERATIONS> sequences(NUM_RANDOM_ITERATIONS, "");
  for (usize idx = 0; idx < NUM_RANDOM_ITERATIONS; ++idx) {
    sequences.at(idx) = GenerateRandomDnaSequence(SEQ_LEN);
  }

  for (auto const& sequence : sequences) {
    auto const slides_list = SlidingView(sequence, 12);
    for (auto const& seq : slides_list) {
      auto const rc_seq = RevComp(seq);

      CAPTURE(seq.substr(0, KMER_SIZE), seq.substr(1, KMER_SIZE));
      CAPTURE(rc_seq.substr(0, KMER_SIZE), rc_seq.substr(1, KMER_SIZE));

      SECTION("Forward direction merge") {
        auto first = Kmer(seq.substr(0, KMER_SIZE));
        auto const second = Kmer(seq.substr(1, KMER_SIZE));
        auto const fwd_edge = MakeFwdEdgeKind({first.SignFor(DFLT_ORD), second.SignFor(DFLT_ORD)});

        CAPTURE(first.SequenceFor(DFLT_ORD), first.SignFor(DFLT_ORD));
        CAPTURE(second.SequenceFor(DFLT_ORD), second.SignFor(DFLT_ORD));
        CAPTURE(fwd_edge);

        first.Merge(second, fwd_edge, KMER_SIZE);
        CAPTURE(first.SequenceFor(DFLT_ORD), first.SignFor(DFLT_ORD));
        REQUIRE(MatchesOneOfTwo(first.SequenceFor(DFLT_ORD), {seq, rc_seq}));
      }

      SECTION("Reverse direction merge") {
        auto first = Kmer(seq.substr(1, KMER_SIZE));
        auto const second = Kmer(seq.substr(0, KMER_SIZE));
        auto const rev_edge =
            RevEdgeKind(MakeFwdEdgeKind({second.SignFor(DFLT_ORD), first.SignFor(DFLT_ORD)}));

        CAPTURE(first.SequenceFor(DFLT_ORD), first.SignFor(DFLT_ORD));
        CAPTURE(second.SequenceFor(DFLT_ORD), second.SignFor(DFLT_ORD));
        CAPTURE(rev_edge);

        first.Merge(second, rev_edge, KMER_SIZE);
        CAPTURE(first.SequenceFor(DFLT_ORD), first.SignFor(DFLT_ORD));
        REQUIRE(MatchesOneOfTwo(first.SequenceFor(DFLT_ORD), {seq, rc_seq}));
      }
    }
  }
}

TEST_CASE("Can merge two adjacent unequal sized kmers", "[lancet][cbdg][Kmer]") {
  std::random_device device;
  std::mt19937_64 generator(device());

  static constexpr usize MIN_KMER_SIZE = 11;
  static constexpr usize MAX_KMER_SIZE = 101;
  static constexpr usize MAX_SEQ_LEN = 999;

  std::uniform_int_distribution<usize> seq_len_picker(3 * MAX_KMER_SIZE, MAX_SEQ_LEN);
  std::uniform_int_distribution<usize> ksize_picker(MIN_KMER_SIZE, MAX_KMER_SIZE);

  for (usize idx = 0; idx < NUM_RANDOM_ITERATIONS; ++idx) {
    auto const kmer_size = static_cast<usize>((2 * ksize_picker(generator)) + 1);
    auto const total_length = static_cast<usize>(seq_len_picker(generator));

    auto const first_length = kmer_size;
    auto const second_start = first_length - kmer_size + 1;
    auto const second_length = (total_length - first_length) + (kmer_size - 1);

    CAPTURE(kmer_size, total_length, first_length, second_start, second_length);
    auto const sequence = GenerateRandomDnaSequence(total_length);
    std::string_view const seq = sequence;
    auto const rc_seq = RevComp(sequence);

    SECTION("Forward direction merge") {
      auto fwd_first = Kmer(seq.substr(0, first_length));
      auto const fwd_second = Kmer(seq.substr(second_start, second_length));
      auto const fwd_edge =
          MakeFwdEdgeKind({fwd_first.SignFor(DFLT_ORD), fwd_second.SignFor(DFLT_ORD)});

      CAPTURE(fwd_first.SequenceFor(DFLT_ORD), fwd_first.SignFor(DFLT_ORD));
      CAPTURE(fwd_second.SequenceFor(DFLT_ORD), fwd_second.SignFor(DFLT_ORD));
      CAPTURE(fwd_edge);

      fwd_first.Merge(fwd_second, fwd_edge, kmer_size);
      CAPTURE(fwd_first.SequenceFor(DFLT_ORD), fwd_first.SignFor(DFLT_ORD));
      REQUIRE(MatchesOneOfTwo(fwd_first.SequenceFor(DFLT_ORD), {seq, rc_seq}));
    }

    SECTION("Reverse direction merge") {
      auto rev_first = Kmer(seq.substr(second_start, second_length));
      auto const rev_second = Kmer(seq.substr(0, first_length));
      auto const rev_edge =
          RevEdgeKind(MakeFwdEdgeKind({rev_second.SignFor(DFLT_ORD), rev_first.SignFor(DFLT_ORD)}));

      CAPTURE(rev_first.SequenceFor(DFLT_ORD), rev_first.SignFor(DFLT_ORD));
      CAPTURE(rev_second.SequenceFor(DFLT_ORD), rev_second.SignFor(DFLT_ORD));
      CAPTURE(rev_edge);

      rev_first.Merge(rev_second, rev_edge, kmer_size);
      CAPTURE(rev_first.SequenceFor(DFLT_ORD), rev_first.SignFor(DFLT_ORD));
      REQUIRE(MatchesOneOfTwo(rev_first.SequenceFor(DFLT_ORD), {seq, rc_seq}));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Can merge multiple adjacent equal sized kmers", "[lancet][cbdg][Kmer]") {
  for (usize idx = 0; idx < NUM_RANDOM_ITERATIONS; ++idx) {
    static constexpr usize LONG_SEQ_LEN = 1024;
    static constexpr usize MER_SIZE = 21;

    auto const sequence = GenerateRandomDnaSequence(LONG_SEQ_LEN);
    auto const rc_sequence = RevComp(sequence);
    auto const mers_list = SlidingKmers(sequence, MER_SIZE);

    SECTION("Forward direction merge") {
      Kmer merged_seq;
      for (auto const& mer : mers_list) {
        auto const fwd_edge =
            MakeFwdEdgeKind({merged_seq.SignFor(DFLT_ORD), mer.SignFor(DFLT_ORD)});
        CAPTURE(merged_seq.SequenceFor(DFLT_ORD), merged_seq.SignFor(DFLT_ORD));
        CAPTURE(mer.SequenceFor(DFLT_ORD), mer.SignFor(DFLT_ORD));
        CAPTURE(fwd_edge);
        merged_seq.Merge(mer, fwd_edge, MER_SIZE);
      }

      CHECK(MatchesOneOfTwo(merged_seq.SequenceFor(DFLT_ORD), {sequence, rc_sequence}));
    }

    SECTION("Reverse direction merge") {
      Kmer rev_merged_seq;
      for (auto const& mer : std::ranges::reverse_view(mers_list)) {
        auto const rev_edge =
            RevEdgeKind(MakeFwdEdgeKind({mer.SignFor(DFLT_ORD), rev_merged_seq.SignFor(DFLT_ORD)}));
        CAPTURE(rev_merged_seq.SequenceFor(DFLT_ORD), rev_merged_seq.SignFor(DFLT_ORD));
        CAPTURE(mer.SequenceFor(DFLT_ORD), mer.SignFor(DFLT_ORD));
        CAPTURE(rev_edge);
        rev_merged_seq.Merge(mer, rev_edge, MER_SIZE);
      }

      CHECK(MatchesOneOfTwo(rev_merged_seq.SequenceFor(DFLT_ORD), {sequence, rc_sequence}));
    }
  }
}
