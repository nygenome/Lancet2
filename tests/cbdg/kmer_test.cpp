#include "lancet/cbdg/kmer.h"

#include "lancet/base/assert.h"
#include "lancet/base/rev_comp.h"
#include "lancet/base/sliding.h"
#include "lancet/base/types.h"

#include "absl/container/fixed_array.h"
#include "absl/random/distributions.h"
#include "absl/strings/string_view.h"
#include "catch_amalgamated.hpp"

#include <array>
#include <iterator>
#include <random>
#include <ranges>
#include <string>
#include <string_view>

using lancet::cbdg::Kmer;
using lancet::cbdg::MakeFwdEdgeKind;
using lancet::cbdg::RevEdgeKind;

namespace {

/// Generate a random DNA sequence of the given length for fuzz testing.
/// Caller passes a const literal seed so the same TEST_CASE produces
/// the same sequence on every run — `std::random_device` is forbidden
/// in tests (per the project's determinism convention) because two
/// failing runs would be impossible to compare.
inline auto GenerateRandomDnaSequence(usize const seq_len, u64 const seed) -> std::string {
  static constexpr std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};

  // Const-literal seed is the project's documented determinism
  // convention. The clang-tidy check is conservatively designed for
  // production code; in tests, predictability is exactly what we want.
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937_64 generator(seed);

  std::string result(seq_len, 'N');

  for (usize iter = 0; iter < seq_len; ++iter) {
    result[iter] = BASES.at(absl::Uniform<usize>(absl::IntervalClosed, generator, 0, 3));
  }

  return result;
}

/// Return true if `result` matches either of the two expected values (fwd or revcomp).
inline auto MatchesOneOfTwo(std::string_view result, std::array<std::string_view, 2> const& values)
    -> bool {
  return (result == values[0]) || (result == values[1]);
}

/// Decompose a sequence into a sliding window of Kmers for merge testing.
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
// Pinned base seed: every iteration derives a distinct seed by adding
// the loop index, so iterations produce deterministic sequences across
// runs. Re-running a TEST_CASE on a different machine (or after an
// unrelated refactor) reproduces the same sequences on every iteration.
static constexpr u64 BASE_SEED = 0x5E'ED'5E'ED'5E'ED'5E'EDULL;

// Catch2 SECTION fan-out inflates clang-tidy's cognitive-complexity metric beyond the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Can merge two adjacent equal sized kmers", "[lancet][cbdg][Kmer]") {
  static constexpr auto SEQ_LEN = 151;
  static constexpr auto KMER_SIZE = 11;

  absl::FixedArray<std::string, NUM_RANDOM_ITERATIONS> sequences(NUM_RANDOM_ITERATIONS, "");
  for (usize iter = 0; iter < NUM_RANDOM_ITERATIONS; ++iter) {
    sequences.at(iter) = GenerateRandomDnaSequence(SEQ_LEN, BASE_SEED + iter);
  }

  for (auto const& sequence : sequences) {
    auto const slides_list = lancet::base::SlidingView(sequence, 12);
    for (auto const& slide : slides_list) {
      auto const rc_seq = lancet::base::RevComp(slide);

      CAPTURE(slide.substr(0, KMER_SIZE), slide.substr(1, KMER_SIZE));
      CAPTURE(rc_seq.substr(0, KMER_SIZE), rc_seq.substr(1, KMER_SIZE));

      SECTION("Forward direction merge") {
        auto first = Kmer(slide.substr(0, KMER_SIZE));
        auto const second = Kmer(slide.substr(1, KMER_SIZE));
        auto const fwd_edge = MakeFwdEdgeKind({first.SignFor(DFLT_ORD), second.SignFor(DFLT_ORD)});

        CAPTURE(first.SequenceFor(DFLT_ORD), first.SignFor(DFLT_ORD));
        CAPTURE(second.SequenceFor(DFLT_ORD), second.SignFor(DFLT_ORD));
        CAPTURE(fwd_edge);

        first.Merge(second, fwd_edge, KMER_SIZE);
        CAPTURE(first.SequenceFor(DFLT_ORD), first.SignFor(DFLT_ORD));
        REQUIRE(MatchesOneOfTwo(first.SequenceFor(DFLT_ORD), {slide, rc_seq}));
      }

      SECTION("Reverse direction merge") {
        auto first = Kmer(slide.substr(1, KMER_SIZE));
        auto const second = Kmer(slide.substr(0, KMER_SIZE));
        auto const rev_edge =
            RevEdgeKind(MakeFwdEdgeKind({second.SignFor(DFLT_ORD), first.SignFor(DFLT_ORD)}));

        CAPTURE(first.SequenceFor(DFLT_ORD), first.SignFor(DFLT_ORD));
        CAPTURE(second.SequenceFor(DFLT_ORD), second.SignFor(DFLT_ORD));
        CAPTURE(rev_edge);

        first.Merge(second, rev_edge, KMER_SIZE);
        CAPTURE(first.SequenceFor(DFLT_ORD), first.SignFor(DFLT_ORD));
        REQUIRE(MatchesOneOfTwo(first.SequenceFor(DFLT_ORD), {slide, rc_seq}));
      }
    }
  }
}

TEST_CASE("Can merge two adjacent unequal sized kmers", "[lancet][cbdg][Kmer]") {
  // Const-literal seed is the project's documented determinism
  // convention. The clang-tidy check is conservatively designed for
  // production code; in tests, predictability is exactly what we want.
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937_64 generator(BASE_SEED);

  static constexpr usize MIN_KMER_SIZE = 11;
  static constexpr usize MAX_KMER_SIZE = 101;
  static constexpr usize MAX_SEQ_LEN = 999;

  for (usize iter = 0; iter < NUM_RANDOM_ITERATIONS; ++iter) {
    auto const kmer_size =
        (2 * absl::Uniform<usize>(absl::IntervalClosed, generator, MIN_KMER_SIZE, MAX_KMER_SIZE)) +
        1;
    auto const total_length =
        absl::Uniform<usize>(absl::IntervalClosed, generator, 3 * MAX_KMER_SIZE, MAX_SEQ_LEN);

    auto const first_length = kmer_size;
    auto const second_start = first_length - kmer_size + 1;
    auto const second_length = (total_length - first_length) + (kmer_size - 1);

    CAPTURE(kmer_size, total_length, first_length, second_start, second_length);
    auto const sequence = GenerateRandomDnaSequence(total_length, BASE_SEED + iter);
    std::string_view const slide = sequence;
    auto const rc_seq = lancet::base::RevComp(sequence);

    SECTION("Forward direction merge") {
      auto fwd_first = Kmer(slide.substr(0, first_length));
      auto const fwd_second = Kmer(slide.substr(second_start, second_length));
      auto const fwd_edge =
          MakeFwdEdgeKind({fwd_first.SignFor(DFLT_ORD), fwd_second.SignFor(DFLT_ORD)});

      CAPTURE(fwd_first.SequenceFor(DFLT_ORD), fwd_first.SignFor(DFLT_ORD));
      CAPTURE(fwd_second.SequenceFor(DFLT_ORD), fwd_second.SignFor(DFLT_ORD));
      CAPTURE(fwd_edge);

      fwd_first.Merge(fwd_second, fwd_edge, kmer_size);
      CAPTURE(fwd_first.SequenceFor(DFLT_ORD), fwd_first.SignFor(DFLT_ORD));
      REQUIRE(MatchesOneOfTwo(fwd_first.SequenceFor(DFLT_ORD), {slide, rc_seq}));
    }

    SECTION("Reverse direction merge") {
      auto rev_first = Kmer(slide.substr(second_start, second_length));
      auto const rev_second = Kmer(slide.substr(0, first_length));
      auto const rev_edge =
          RevEdgeKind(MakeFwdEdgeKind({rev_second.SignFor(DFLT_ORD), rev_first.SignFor(DFLT_ORD)}));

      CAPTURE(rev_first.SequenceFor(DFLT_ORD), rev_first.SignFor(DFLT_ORD));
      CAPTURE(rev_second.SequenceFor(DFLT_ORD), rev_second.SignFor(DFLT_ORD));
      CAPTURE(rev_edge);

      rev_first.Merge(rev_second, rev_edge, kmer_size);
      CAPTURE(rev_first.SequenceFor(DFLT_ORD), rev_first.SignFor(DFLT_ORD));
      REQUIRE(MatchesOneOfTwo(rev_first.SequenceFor(DFLT_ORD), {slide, rc_seq}));
    }
  }
}

// Catch2 SECTION fan-out inflates clang-tidy's cognitive-complexity metric beyond the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Can merge multiple adjacent equal sized kmers", "[lancet][cbdg][Kmer]") {
  for (usize iter = 0; iter < NUM_RANDOM_ITERATIONS; ++iter) {
    static constexpr usize LONG_SEQ_LEN = 1024;
    static constexpr usize MER_SIZE = 21;

    auto const sequence = GenerateRandomDnaSequence(LONG_SEQ_LEN, BASE_SEED + iter);
    auto const rc_sequence = lancet::base::RevComp(sequence);
    auto const mers_list = SlidingKmers(sequence, MER_SIZE);

    SECTION("Forward direction merge") {
      Kmer merged_seq;
      for (auto const& kmer : mers_list) {
        auto const fwd_edge =
            MakeFwdEdgeKind({merged_seq.SignFor(DFLT_ORD), kmer.SignFor(DFLT_ORD)});
        CAPTURE(merged_seq.SequenceFor(DFLT_ORD), merged_seq.SignFor(DFLT_ORD));
        CAPTURE(kmer.SequenceFor(DFLT_ORD), kmer.SignFor(DFLT_ORD));
        CAPTURE(fwd_edge);
        merged_seq.Merge(kmer, fwd_edge, MER_SIZE);
      }

      CHECK(MatchesOneOfTwo(merged_seq.SequenceFor(DFLT_ORD), {sequence, rc_sequence}));
    }

    SECTION("Reverse direction merge") {
      Kmer rev_merged_seq;
      for (auto const& kmer : std::ranges::reverse_view(mers_list)) {
        auto const rev_edge = RevEdgeKind(
            MakeFwdEdgeKind({kmer.SignFor(DFLT_ORD), rev_merged_seq.SignFor(DFLT_ORD)}));
        CAPTURE(rev_merged_seq.SequenceFor(DFLT_ORD), rev_merged_seq.SignFor(DFLT_ORD));
        CAPTURE(kmer.SequenceFor(DFLT_ORD), kmer.SignFor(DFLT_ORD));
        CAPTURE(rev_edge);
        rev_merged_seq.Merge(kmer, rev_edge, MER_SIZE);
      }

      CHECK(MatchesOneOfTwo(rev_merged_seq.SequenceFor(DFLT_ORD), {sequence, rc_sequence}));
    }
  }
}
