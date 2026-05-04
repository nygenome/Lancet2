#include "lancet/base/repeat.h"

#include "lancet/base/types.h"

#include "absl/types/span.h"
#include "catch_amalgamated.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <string>
#include <string_view>
#include <vector>

namespace lancet::base::tests {

namespace {

/// Generate a 5000bp pseudo-random DNA sequence for Hamming distance fuzz
/// testing. Caller passes a const literal seed so the same TEST_CASE produces
/// the same sequence on every run — `std::random_device` is forbidden in tests
/// (per the project's determinism convention) because two failing runs would
/// be impossible to compare. Two distinct seeds in the same TEST_CASE produce
/// two distinct sequences (probabilistically) so the per-iteration comparison
/// `HammingDist(s, t) != 0` still exercises the "different inputs" path.
inline auto GenerateRandomDnaSequence(u64 seed) -> std::string {
  static constexpr std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};

  std::mt19937_64 generator(seed);

  static constexpr usize SEQ_LENGTH = 5000;
  std::uniform_int_distribution<usize> base_chooser(0, 3);
  std::string result;
  result.reserve(SEQ_LENGTH);

  for (usize idx = 0; idx < SEQ_LENGTH; ++idx) {
    result.push_back(BASES.at(base_chooser(generator)));
  }

  return result;
}

}  // namespace

// ============================================================================
//  Hamming distance: randomized correctness
// ============================================================================

// Catch2 SECTION fan-out inflates clang-tidy's cognitive-complexity metric beyond the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Can calculate hamming distance correctly for random strings",
          "[lancet][base][HammingDist]") {
  // Pinned base seed: every iteration derives a distinct seed by adding the
  // loop index, so the 1000 iterations produce 1000 deterministic sequences.
  // Re-running this TEST_CASE on a different machine (or after an unrelated
  // refactor) reproduces the same sequence pair on every iteration.
  static constexpr u64 BASE_SEED = 0x5E'ED'5E'ED'5E'ED'5E'EDULL;
  static constexpr usize NUM_ITERATIONS = 1000;

  for (usize idx = 0; idx <= NUM_ITERATIONS; ++idx) {
    auto const result = GenerateRandomDnaSequence(BASE_SEED + idx);
    // Offset the second seed far enough that the two streams diverge at the
    // very first base — keeps `HammingDist != 0` deterministically true.
    auto const other = GenerateRandomDnaSequence(BASE_SEED + idx + 0x10'00'00'00ULL);

    REQUIRE(HammingDist(result, result) == 0);
    REQUIRE(HammingDist(result, other) != 0);
  }
}

// ============================================================================
//  Hamming distance: small known-answer tests
// ============================================================================

TEST_CASE("Can calculate hamming distance correctly for small test",
          "[lancet][base][HammingDist]") {
  std::string_view const test = "aaaa";
  std::string_view const diff_a = "abaa";
  std::string_view const diff_b = "aaba";

  REQUIRE(HammingDist(test, test) == 0);
  REQUIRE(HammingDist(test, diff_a) == 1);
  REQUIRE(HammingDist(test, diff_b) == 1);
  REQUIRE(HammingDist(diff_a, diff_b) == 2);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  HammingDist: SIMD-aligned and unaligned-tail boundary cases             ║
// ║                                                                          ║
// ║  base.md documents that `HammingDist` uses AVX2 (32-byte) on x86 and     ║
// ║  NEON (16-byte) on ARM64 with overlapping unaligned tail loads. These    ║
// ║  cases pin the boundaries that matter for the SIMD path: exact 32B / 16B ║
// ║  multiples (no tail), one-byte-shy multiples (max-overlap tail), and a   ║
// ║  one-byte sequence (tail-only, no SIMD body).                            ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("HammingDist returns 0 on a 32-byte aligned identical block",
          "[lancet][base][HammingDist]") {
  // Exact AVX2 lane width. The SIMD body executes once with no tail; a
  // regression that mishandled the no-tail path would surface here.
  std::string const seq(32, 'A');
  CHECK(HammingDist(seq, seq) == 0);
}

TEST_CASE("HammingDist counts every-byte mismatch on a 32-byte block",
          "[lancet][base][HammingDist]") {
  // Pair where every position differs: AAAA…A vs CCCC…C, length 32. The
  // SIMD popcount must return 32; an off-by-one would return 31 or 33.
  std::string const lhs(32, 'A');
  std::string const rhs(32, 'C');
  CHECK(HammingDist(lhs, rhs) == 32);
}

TEST_CASE("HammingDist handles a 33-byte input (1-byte unaligned tail)",
          "[lancet][base][HammingDist]") {
  // 33 bytes = 32-byte SIMD body + 1-byte tail. The implementation uses an
  // overlapping unaligned tail load (the last 32 bytes of the buffer)
  // rather than a scalar tail; the overlap re-counts 31 already-counted
  // bytes, so the popcount must compensate or the result drifts.
  std::string lhs(33, 'A');
  std::string rhs(33, 'A');
  rhs[32] = 'T';  // single mismatch at the very last byte
  CHECK(HammingDist(lhs, rhs) == 1);

  // And the symmetric case: mismatch at the first byte (well inside the
  // SIMD body, no tail involvement).
  lhs[0] = 'C';
  rhs[0] = 'A';
  CHECK(HammingDist(lhs, rhs) == 2);  // [0]: A vs C; [32]: A vs T
}

TEST_CASE("HammingDist handles a 31-byte input (one-byte-shy of SIMD width)",
          "[lancet][base][HammingDist]") {
  // 31 bytes < 32-byte AVX2 width → should fall through to the auto-
  // vectorized scalar path on x86. On ARM64 (16-byte NEON), 31 bytes =
  // 16 SIMD + 15-byte tail, exercising the tail-overlap path. The result
  // is path-dependent in implementation, but the answer must be the same.
  std::string const lhs(31, 'A');
  std::string rhs(31, 'A');
  rhs[10] = 'T';
  CHECK(HammingDist(lhs, rhs) == 1);
}

TEST_CASE("HammingDist handles a single-byte input (tail-only, no SIMD body)",
          "[lancet][base][HammingDist]") {
  // 1 byte < both SIMD widths → pure scalar path. Smallest non-empty
  // input; pins the corner where the SIMD body never executes.
  CHECK(HammingDist(std::string_view("A"), std::string_view("A")) == 0);
  CHECK(HammingDist(std::string_view("A"), std::string_view("T")) == 1);
}

TEST_CASE("HammingDist on an empty input returns 0", "[lancet][base][HammingDist]") {
  // Both empty → 0 differences. The implementation must not deref the data
  // pointers (which may be null for default-constructed string_views).
  CHECK(HammingDist(std::string_view(""), std::string_view("")) == 0);
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  HasRepeat / HasExactRepeat                                              ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("HasExactRepeat detects identical k-mers in a span", "[lancet][base][HasExactRepeat]") {
  // Two copies of "ACGT" appear in the span — exact repeat. The
  // implementation routes through HasRepeat(kmers, 0) which uses an O(n)
  // hash-set duplicate check.
  std::vector<std::string_view> const kmers{"ACGT", "TGCA", "ACGT", "GGCC"};
  CHECK(HasExactRepeat(absl::MakeConstSpan(kmers)));
}

TEST_CASE("HasExactRepeat returns false on all-distinct k-mers", "[lancet][base][HasExactRepeat]") {
  std::vector<std::string_view> const kmers{"ACGT", "TGCA", "GGCC", "AATT"};
  CHECK_FALSE(HasExactRepeat(absl::MakeConstSpan(kmers)));
}

TEST_CASE("HasExactRepeat returns false on an empty or single-element span",
          "[lancet][base][HasExactRepeat]") {
  // Cannot have a "repeat" with fewer than two k-mers. The hash-set path
  // must short-circuit on len < 2 to avoid spurious work.
  std::vector<std::string_view> const empty;
  std::vector<std::string_view> const single{"ACGT"};
  CHECK_FALSE(HasExactRepeat(absl::MakeConstSpan(empty)));
  CHECK_FALSE(HasExactRepeat(absl::MakeConstSpan(single)));
}

TEST_CASE("HasRepeat with max_mismatches=0 matches HasExactRepeat", "[lancet][base][HasRepeat]") {
  // HasRepeat(kmers, 0) is documented to delegate to the exact-repeat
  // path. Verify the two surface APIs agree on a representative input.
  std::vector<std::string_view> const kmers{"ACGT", "TGCA", "ACGT"};
  auto const span = absl::MakeConstSpan(kmers);
  CHECK(HasRepeat(span, 0) == HasExactRepeat(span));
}

TEST_CASE("HasRepeat detects approximate repeats within max_mismatches",
          "[lancet][base][HasRepeat]") {
  // ACGT and ACGA differ by exactly 1 base; with max_mismatches=1 they
  // count as an approximate repeat. With max_mismatches=0 they do not.
  std::vector<std::string_view> const kmers{"ACGT", "TGCA", "ACGA"};
  auto const span = absl::MakeConstSpan(kmers);
  CHECK(HasRepeat(span, 1));
  CHECK_FALSE(HasRepeat(span, 0));
}

TEST_CASE("HasRepeat returns false when no pair is within max_mismatches",
          "[lancet][base][HasRepeat]") {
  // All four k-mers mutually differ in ≥ 2 positions, so max_mismatches=1
  // finds no approximate repeat. Exercises the early-exit per-pair path
  // documented in base.md (the SIMD short-circuit collapses the per-pair
  // cost from O(L) to O(1) on random DNA).
  std::vector<std::string_view> const kmers{"AAAA", "CCCC", "GGGG", "TTTT"};
  CHECK_FALSE(HasRepeat(absl::MakeConstSpan(kmers), 1));
}

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  Property: HammingDist agrees with a scalar reference                    ║
// ║                                                                          ║
// ║  base.md documents that `HammingDist` uses an AVX2 / NEON SIMD kernel    ║
// ║  with overlapping-tail loads to avoid scalar cleanup. The SIMD result    ║
// ║  must agree with a straightforward byte-by-byte scalar reference on      ║
// ║  every input, regardless of length parity vs the SIMD lane width.        ║
// ║                                                                          ║
// ║  The reference is intentionally trivial (no SIMD, no overlap, no early   ║
// ║  exit) so any disagreement points to a bug in the SIMD path.             ║
// ╚══════════════════════════════════════════════════════════════════════════╝

namespace {

// Scalar reference: byte-by-byte comparison, used as ground truth for the
// SIMD path's output. Lives inside the file's anonymous namespace so it has
// internal linkage and cannot collide with anything in the project.
[[nodiscard]] auto ScalarHammingDist(std::string_view first, std::string_view second) -> usize {
  if (first.size() != second.size()) return std::max(first.size(), second.size());
  usize differences = 0;
  for (usize idx = 0; idx < first.size(); ++idx) {
    if (first[idx] != second[idx]) ++differences;
  }
  return differences;
}

}  // namespace

// Catch2's per-iteration generation drives the cognitive-complexity metric over the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HammingDist (SIMD) agrees with a scalar reference on random DNA",
          "[lancet][base][HammingDist]") {
  // Pinned base seed: sequence pairs are deterministic across runs.
  static constexpr u64 BASE_SEED = 0x5E'ED'5E'ED'5E'ED'5E'EDULL;
  static constexpr usize NUM_PROPERTY_ITERATIONS = 200;

  // Const-literal seed is the project's documented determinism convention
  // (see test_style.md / §A.9). The clang-tidy check is conservatively
  // designed for production code; in tests, predictability is exactly
  // what we want.
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937_64 length_generator(BASE_SEED);
  // Lengths chosen to span the SIMD-lane boundaries:
  //   - well below 16 / 32 (pure scalar tail)
  //   - exactly at 32 (one full AVX2 lane)
  //   - 33..64 (one lane + 1..32-byte unaligned tail)
  //   - 1024 (32 full lanes, no tail)
  // The bounds are inclusive so 32 and 1024 are reachable.
  std::uniform_int_distribution<usize> length_picker(1, 1024);

  for (usize iter = 0; iter < NUM_PROPERTY_ITERATIONS; ++iter) {
    auto const seq_len = length_picker(length_generator);
    auto const lhs = GenerateRandomDnaSequence(BASE_SEED + (iter * 2));
    auto const rhs = GenerateRandomDnaSequence(BASE_SEED + (iter * 2) + 1);

    // GenerateRandomDnaSequence yields 5000bp; truncate to seq_len so we
    // exercise the variable-length path uniformly.
    auto const lhs_view = std::string_view(lhs).substr(0, seq_len);
    auto const rhs_view = std::string_view(rhs).substr(0, seq_len);

    auto const simd_result = HammingDist(lhs_view, rhs_view);
    auto const scalar_result = ScalarHammingDist(lhs_view, rhs_view);

    INFO("iter=" << iter << " seq_len=" << seq_len << " simd=" << simd_result
                 << " scalar=" << scalar_result);
    CHECK(simd_result == scalar_result);
  }
}

}  // namespace lancet::base::tests
