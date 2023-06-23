#include "lancet/base/hash.h"

#include "wyhash.h"

// NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
auto HashStr(std::string_view str) -> u64 { return wyhash(str.data(), str.length(), 0, _wyp); }
