#include "lancet/base/hash.h"

#include "MurmurHash3.h"
#include "wyhash.h"

auto HashStr64(std::string_view str) -> u64 {
  static constexpr u64 PRIME_SEED = 18446744073709551557ULL;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
  return wyhash(str.data(), str.length(), PRIME_SEED, _wyp);
}

auto HashStr32(std::string_view str) -> u32 {
  static constexpr u32 PRIME_SEED = 2147483647;
  u32 hash = 0;
  MurmurHash3_x86_32(str.data(), static_cast<int>(str.length()), PRIME_SEED, &hash);
  return hash;
}
