#include "lancet/base/repeat.h"

#include "lancet/base/assert.h"
#include "lancet/base/types.h"

#include "absl/container/flat_hash_set.h"
#include "absl/types/span.h"

#include <algorithm>
#include <string_view>
#include <utility>

// Platform-specific SIMD intrinsics — CMake guarantees -march=x86-64-v3 on x86
// (AVX2 + BMI2 + POPCNT) and aarch64 baseline NEON on ARM64.
#ifdef __AVX2__
#include <emmintrin.h>
#include <immintrin.h>
#include <mmintrin.h>
#elif defined(__aarch64__) || defined(_M_ARM64)
#include <arm_neon.h>
#endif

namespace lancet::base {

namespace {

// ============================================================================
// IsWithinHammingDist — SIMD-accelerated early-exit mismatch check
//
// Returns true if two equal-length byte strings differ in at most
// `max_mismatches` positions.  The key optimization is the early-exit
// after each SIMD chunk: most random DNA pairs exceed a small threshold
// within the first 16–32 bytes, collapsing the per-pair cost from O(L)
// to effectively O(1) in the O(n²) repeat detection loop.
//
// x86 AVX2 path:
//   _mm256_cmpeq_epi8  → byte-wise equality (0xFF match, 0x00 mismatch)
//   _mm256_movemask_epi8 → extract MSB per byte into 32-bit scalar
//   ~mask + __builtin_popcount → count mismatches via hardware POPCNT
//
// ARM64 NEON path:
//   vceqq_u8   → byte-wise equality
//   vbicq_u8   → single-cycle Bit Clear: (ones AND NOT eq) maps
//                 matches → 0, mismatches → 1 without extra inversion
//   vaddvq_u8  → horizontal byte sum across 16 lanes
//
// Tail handling uses overlapping unaligned loads rather than a scalar
// cleanup loop.  On x86 the overlap region is excluded by right-shifting
// the scalar bitmask.  On ARM64 the overlap is zeroed via a precomputed
// TAIL_MASK lookup table indexed by the number of valid bytes.
// ============================================================================

// SIMD paths are flat, not deeply nested
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
[[nodiscard]] inline auto IsWithinHammingDist(std::string_view first, std::string_view second,
                                              usize max_mismatches) -> bool {
  LANCET_ASSERT(first.length() == second.length())

  auto const length = first.length();
  // SIMD load intrinsics require u8*
  // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* ptr1 = reinterpret_cast<u8 const*>(first.data());
  auto const* ptr2 = reinterpret_cast<u8 const*>(second.data());
  // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

  usize dist = 0;
  usize idx = 0;

#ifdef __AVX2__
  // ============================================================================
  // AVX2: 32-byte chunks with early exit
  // ============================================================================
  if (length >= 32) {
    for (; idx + 32 <= length; idx += 32) {
      // SIMD load intrinsic _mm256_loadu_si256 requires `__m256i const*` per Intel API contract.
      // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const vec1 = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(ptr1 + idx));
      auto const vec2 = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(ptr2 + idx));
      // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const cmp = _mm256_cmpeq_epi8(vec1, vec2);

      // movemask: MSB of each byte → 32-bit scalar.  Equal bytes have MSB=1, so
      // ~mask flips to mismatches=1.  __builtin_popcount compiles to hardware POPCNT
      // (guaranteed by -march=x86-64-v3 which implies BMI2 + POPCNT).
      auto const mask = ~static_cast<u32>(_mm256_movemask_epi8(cmp));
      dist += static_cast<usize>(__builtin_popcount(mask));
      if (dist > max_mismatches) return false;
    }

    // Overlapping tail: re-process the last 32 bytes, masking out already-counted
    // overlap via right-shift.  Avoids a scalar cleanup loop entirely.
    if (idx < length) {
      auto const offset = length - 32;
      // SIMD load intrinsic _mm256_loadu_si256 requires `__m256i const*` per Intel API contract.
      // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const vec1 = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(ptr1 + offset));
      auto const vec2 = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(ptr2 + offset));
      // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const cmp = _mm256_cmpeq_epi8(vec1, vec2);
      auto mask = ~static_cast<u32>(_mm256_movemask_epi8(cmp));

      // Lower `overlap` bits correspond to bytes already counted in the main loop.
      // Right-shift discards them, keeping only the unprocessed tail bits.
      auto const overlap = static_cast<u32>(32 - (length - idx));
      mask >>= overlap;
      dist += static_cast<usize>(__builtin_popcount(mask));
    }
    return dist <= max_mismatches;
  }

  // ============================================================================
  // SSE: 16-byte path for lengths 16–31
  // ============================================================================
  if (length >= 16) {
    // SIMD load intrinsic _mm_loadu_si128 requires `__m128i const*` per Intel API contract.
    // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
    auto const vec1 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(ptr1));
    auto const vec2 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(ptr2));
    // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
    auto const cmp = _mm_cmpeq_epi8(vec1, vec2);

    // _mm_movemask_epi8 returns 16-bit result in 32-bit int.  The & 0xFFFF
    // prevents ~mask from setting upper 16 bits to 1, which would inflate popcount.
    auto mask = ~static_cast<u32>(_mm_movemask_epi8(cmp)) & 0xFFFF;
    dist += static_cast<usize>(__builtin_popcount(mask));
    if (dist > max_mismatches) return false;

    idx = 16;
    if (idx < length) {
      auto const offset = length - 16;
      // SIMD load intrinsic _mm_loadu_si128 requires `__m128i const*` per Intel API contract.
      // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const vec1_tail = _mm_loadu_si128(reinterpret_cast<__m128i const*>(ptr1 + offset));
      auto const vec2_tail = _mm_loadu_si128(reinterpret_cast<__m128i const*>(ptr2 + offset));
      // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const eq_tail = _mm_cmpeq_epi8(vec1_tail, vec2_tail);
      auto mask_tail = ~static_cast<u32>(_mm_movemask_epi8(eq_tail)) & 0xFFFF;

      auto const overlap = static_cast<u32>(16 - (length - idx));
      mask_tail >>= overlap;
      dist += static_cast<usize>(__builtin_popcount(mask_tail));
    }
    return dist <= max_mismatches;
  }

#elif defined(__aarch64__) || defined(_M_ARM64)
  // ============================================================================
  // NEON: 16-byte chunks with early exit
  // ============================================================================
  if (length >= 16) {
    auto const ones = vdupq_n_u8(1);
    for (; idx + 16 <= length; idx += 16) {
      auto const vec1 = vld1q_u8(ptr1 + idx);
      auto const vec2 = vld1q_u8(ptr2 + idx);
      auto const eq = vceqq_u8(vec1, vec2);

      // vbicq_u8(ones, eq) = ones AND NOT eq:
      //   match   (eq=0xFF) → 0x01 & 0x00 = 0  (no mismatch)
      //   mismatch(eq=0x00) → 0x01 & 0xFF = 1  (mismatch counted)
      // vaddvq_u8 horizontally sums all 16 lanes.  Max per chunk = 16, fits in u8.
      dist += static_cast<usize>(vaddvq_u8(vbicq_u8(ones, eq)));
      if (dist > max_mismatches) return false;
    }

    // Overlapping tail: load from (length - 16), then zero-out already-counted lanes
    // via a precomputed mask indexed by the number of valid (unprocessed) bytes.
    if (idx < length) {
      // clang-format off
      // 32-byte lookup table: TAIL_MASK[valid_bytes] loads 16 bytes where the first
      // (16 - valid_bytes) are 0 (overlap, suppressed) and the last valid_bytes are 1.
      alignas(16) static constexpr u8 TAIL_MASK[32] = {
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
      };
      // clang-format on

      auto const offset = length - 16;
      auto const vec1 = vld1q_u8(ptr1 + offset);
      auto const vec2 = vld1q_u8(ptr2 + offset);
      auto const eq = vceqq_u8(vec1, vec2);
      auto mismatches = vbicq_u8(ones, eq);

      // valid_bytes ranges [1, 15]: TAIL_MASK[valid_bytes] reads bytes
      // [valid_bytes, valid_bytes+15], max index = 30, within the 32-byte array.
      auto const valid_bytes = length - idx;
      auto const tail_mask = vld1q_u8(&TAIL_MASK[valid_bytes]);
      mismatches = vandq_u8(mismatches, tail_mask);

      dist += static_cast<usize>(vaddvq_u8(mismatches));
    }
    return dist <= max_mismatches;
  }
#endif

  // ============================================================================
  // Portable scalar fallback for short sequences or unknown ISA
  // ============================================================================
  for (; idx < length; ++idx) {
    dist += static_cast<usize>(ptr1[idx] != ptr2[idx]);
    if (dist > max_mismatches) return false;
  }

  return true;
}

}  // namespace

// ============================================================================
// HammingDist — SIMD-accelerated exact mismatch count
//
// Same overlapping-tail intrinsic methodology as IsWithinHammingDist but
// without the early-exit threshold check.  Computes the full distance.
//
// This function is the public API for callers that need the exact count.
// In the hot repeat-detection path, IsWithinHammingDist (above) is used
// instead because its early exit collapses the per-pair cost from O(L)
// to O(1) for random DNA.
// ============================================================================
auto HammingDist(std::string_view first, std::string_view second) -> usize {
  LANCET_ASSERT(first.length() == second.length())

  auto const length = first.length();
  // SIMD load intrinsics require u8*
  // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* ptr1 = reinterpret_cast<u8 const*>(first.data());
  auto const* ptr2 = reinterpret_cast<u8 const*>(second.data());
  // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

  usize dist = 0;
  usize idx = 0;

#ifdef __AVX2__
  if (length >= 32) {
    for (; idx + 32 <= length; idx += 32) {
      // SIMD load intrinsic _mm256_loadu_si256 requires `__m256i const*` per Intel API contract.
      // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const vec1 = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(ptr1 + idx));
      auto const vec2 = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(ptr2 + idx));
      // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const cmp = _mm256_cmpeq_epi8(vec1, vec2);
      auto const mask = ~static_cast<u32>(_mm256_movemask_epi8(cmp));
      dist += static_cast<usize>(__builtin_popcount(mask));
    }
    if (idx < length) {
      auto const offset = length - 32;
      // SIMD load intrinsic _mm256_loadu_si256 requires `__m256i const*` per Intel API contract.
      // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const vec1 = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(ptr1 + offset));
      auto const vec2 = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(ptr2 + offset));
      // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const cmp = _mm256_cmpeq_epi8(vec1, vec2);
      auto mask = ~static_cast<u32>(_mm256_movemask_epi8(cmp));
      auto const overlap = static_cast<u32>(32 - (length - idx));
      mask >>= overlap;
      dist += static_cast<usize>(__builtin_popcount(mask));
    }
    return dist;
  }

  if (length >= 16) {
    // SIMD load intrinsic _mm_loadu_si128 requires `__m128i const*` per Intel API contract.
    // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
    auto const vec1 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(ptr1));
    auto const vec2 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(ptr2));
    // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
    auto const cmp = _mm_cmpeq_epi8(vec1, vec2);
    auto mask = ~static_cast<u32>(_mm_movemask_epi8(cmp)) & 0xFFFF;
    dist += static_cast<usize>(__builtin_popcount(mask));

    idx = 16;
    if (idx < length) {
      auto const offset = length - 16;
      // SIMD load intrinsic _mm_loadu_si128 requires `__m128i const*` per Intel API contract.
      // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const vec1_tail = _mm_loadu_si128(reinterpret_cast<__m128i const*>(ptr1 + offset));
      auto const vec2_tail = _mm_loadu_si128(reinterpret_cast<__m128i const*>(ptr2 + offset));
      // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
      auto const eq_tail = _mm_cmpeq_epi8(vec1_tail, vec2_tail);
      auto mask_tail = ~static_cast<u32>(_mm_movemask_epi8(eq_tail)) & 0xFFFF;
      auto const overlap = static_cast<u32>(16 - (length - idx));
      mask_tail >>= overlap;
      dist += static_cast<usize>(__builtin_popcount(mask_tail));
    }
    return dist;
  }

#elif defined(__aarch64__) || defined(_M_ARM64)
  if (length >= 16) {
    auto const ones = vdupq_n_u8(1);
    for (; idx + 16 <= length; idx += 16) {
      auto const vec1 = vld1q_u8(ptr1 + idx);
      auto const vec2 = vld1q_u8(ptr2 + idx);
      auto const eq = vceqq_u8(vec1, vec2);
      dist += static_cast<usize>(vaddvq_u8(vbicq_u8(ones, eq)));
    }

    if (idx < length) {
      // clang-format off
      alignas(16) static constexpr u8 TAIL_MASK[32] = {
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
      };
      // clang-format on
      auto const offset = length - 16;
      auto const vec1 = vld1q_u8(ptr1 + offset);
      auto const vec2 = vld1q_u8(ptr2 + offset);
      auto const eq = vceqq_u8(vec1, vec2);
      auto mismatches = vbicq_u8(ones, eq);

      auto const valid_bytes = length - idx;
      auto const tail_mask = vld1q_u8(&TAIL_MASK[valid_bytes]);
      mismatches = vandq_u8(mismatches, tail_mask);

      dist += static_cast<usize>(vaddvq_u8(mismatches));
    }
    return dist;
  }
#endif

  // Portable fallback: u8 batch accumulation for auto-vectorization.
  // Preserved for architectures where SIMD intrinsics are not available.
  while (idx < length) {
    u8 batch_sum = 0;
    auto const batch_end = std::min(idx + 255, length);
    for (; idx < batch_end; ++idx) {
      batch_sum += static_cast<u8>(ptr1[idx] != ptr2[idx]);
    }
    dist += batch_sum;
  }

  return dist;
}

// ============================================================================
// HasRepeat — repeat detector with fast path for exact matches
//
// For exact repeats (max_mismatches=0), uses an O(n) hash-set duplicate
// check — the set detects identical string_views by content, short-circuiting
// on the first collision.
//
// For approximate repeats, delegates to IsWithinHammingDist which combines
// SIMD-accelerated comparison with early exit.  For random DNA (75% per-base
// mismatch rate), a threshold of max_mismatches=3 is exceeded within ~5 bytes,
// so the first 32-byte AVX2 chunk short-circuits >99% of pairs.  This
// collapses the per-pair cost from O(L) to O(1), making the O(n²) outer
// loop vastly cheaper.
// ============================================================================
auto HasRepeat(absl::Span<std::string_view const> kmers, usize const max_mismatches) -> bool {
  // Exact repeat: O(n) hash-set duplicate detection
  if (max_mismatches == 0) {
    absl::flat_hash_set<std::string_view> seen;
    seen.reserve(kmers.size());
    for (auto const kmer : kmers) {
      if (auto [iter, inserted] = seen.insert(kmer); !inserted) return true;
    }
    return false;
  }

  auto const num_kmers = kmers.size();
  if (num_kmers < 2) return false;

  // Upper-triangle pairwise scan: each pair checked exactly once.
  // IsWithinHammingDist early-exits per chunk, so the typical cost for
  // random DNA at small max_mismatches is ~1 SIMD load per pair.
  for (usize i = 0; i < num_kmers; ++i) {
    for (usize j = i + 1; j < num_kmers; ++j) {
      if (IsWithinHammingDist(kmers[i], kmers[j], max_mismatches)) return true;
    }
  }
  return false;
}

auto HasExactRepeat(absl::Span<std::string_view const> kmers) -> bool {
  return HasRepeat(kmers, 0);
}

}  // namespace lancet::base
