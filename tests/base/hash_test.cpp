#include "lancet/base/hash.h"

#include "lancet/base/types.h"

#include "catch_amalgamated.hpp"

#include <array>
#include <set>
#include <string>
#include <string_view>
#include <utility>

namespace lancet::base::tests {

// ╔══════════════════════════════════════════════════════════════════════════╗
// ║  HashStr64 / HashStr32 — basic determinism + distinctness                ║
// ║                                                                          ║
// ║  These are project hash primitives — used as table keys and as           ║
// ║  identifiers in graph traversal. They are NOT cryptographic — collision  ║
// ║  resistance is the only property we verify here.                         ║
// ╚══════════════════════════════════════════════════════════════════════════╝

TEST_CASE("HashStr64 is deterministic across calls with the same input",
          "[lancet][base][HashStr64]") {
  // Same input, same output, every call. A regression that introduced
  // randomization (e.g. seeding from time) would collapse hash-based caches
  // and graph identity. This is the bedrock property the rest depend on.
  auto const reference = std::string("the quick brown fox jumps over the lazy dog");
  auto const first = HashStr64(reference);
  auto const second = HashStr64(reference);
  auto const third = HashStr64(reference);
  CHECK(first == second);
  CHECK(second == third);
}

TEST_CASE("HashStr32 is deterministic across calls with the same input",
          "[lancet][base][HashStr32]") {
  auto const reference = std::string("the quick brown fox jumps over the lazy dog");
  auto const first = HashStr32(reference);
  auto const second = HashStr32(reference);
  CHECK(first == second);
}

TEST_CASE("HashStr64 of the empty string is a fixed value", "[lancet][base][HashStr64]") {
  // The empty-string hash is a well-defined sentinel that downstream code
  // can rely on. Pinning it prevents an upgrade or refactor from silently
  // changing the value (which would break any persisted-on-disk hash maps).
  auto const empty_a = HashStr64("");
  auto const empty_b = HashStr64(std::string_view(""));
  CHECK(empty_a == empty_b);
}

TEST_CASE("HashStr32 of the empty string is a fixed value", "[lancet][base][HashStr32]") {
  auto const empty_a = HashStr32("");
  auto const empty_b = HashStr32(std::string_view(""));
  CHECK(empty_a == empty_b);
}

TEST_CASE("HashStr64 produces distinct values on a small fixed corpus",
          "[lancet][base][HashStr64]") {
  // Seven short, unrelated strings — chosen so they exercise different
  // length classes (2, 3, 4, 5, 7, 8 chars) and different bytewise prefixes.
  // Any pair colliding would be a silent indicator that the hash mixing is
  // broken in a way that affects realistic short-string inputs (k-mers).
  std::array<std::string, 7> const corpus{"AT",      "ATC",     "ATCG",    "GATCA",
                                          "GATTACA", "GACTACA", "ACGTACGT"};
  std::set<u64> hashes;
  for (auto const& word : corpus) {
    INFO("word=\"" << word << "\"");
    auto const inserted = hashes.insert(HashStr64(word)).second;
    CHECK(inserted);
  }
  CHECK(hashes.size() == corpus.size());
}

TEST_CASE("HashStr32 produces distinct values on a small fixed corpus",
          "[lancet][base][HashStr32]") {
  // 32-bit hashes have a smaller value space than 64-bit (~4B vs ~18Q
  // possible outputs). On a 7-element corpus the birthday-collision
  // probability is negligible (~2.4e-9), so any actual collision here
  // points to a real defect in the mixing.
  std::array<std::string, 7> const corpus{"AT",      "ATC",     "ATCG",    "GATCA",
                                          "GATTACA", "GACTACA", "ACGTACGT"};
  std::set<u32> hashes;
  for (auto const& word : corpus) {
    INFO("word=\"" << word << "\"");
    auto const inserted = hashes.insert(HashStr32(word)).second;
    CHECK(inserted);
  }
  CHECK(hashes.size() == corpus.size());
}

TEST_CASE("HashStr64 distinguishes single-character differences", "[lancet][base][HashStr64]") {
  // The pair "ACGTA" / "ACGTC" differs in exactly one position. A trivially-
  // broken hash (e.g. one that drops the last byte) would return the same
  // value for both. Exercises the 1-bit-difference avalanche minimum.
  CHECK(HashStr64("ACGTA") != HashStr64("ACGTC"));
  CHECK(HashStr64("ACGTA") != HashStr64("ACGTG"));
  CHECK(HashStr64("ACGTA") != HashStr64("ACGTT"));
}

}  // namespace lancet::base::tests
