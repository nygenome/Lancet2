#include "lancet/base/hash.h"

#include <string_view>

#include "absl/hash/internal/city.h"
#include "lancet/base/types.h"

auto HashStr64(std::string_view str) -> u64 {
  // NOLINTNEXTLINE(abseil-no-internal-dependencies)
  return absl::hash_internal::CityHash64(str.data(), str.length());
}

auto HashStr32(std::string_view str) -> u32 {
  // NOLINTNEXTLINE(abseil-no-internal-dependencies)
  return absl::hash_internal::CityHash32(str.data(), str.length());
}
