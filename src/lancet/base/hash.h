#ifndef SRC_LANCET_BASE_HASH_H_
#define SRC_LANCET_BASE_HASH_H_

#include <string_view>

#include "lancet/base/types.h"

[[nodiscard]] auto HashStr64(std::string_view str) -> u64;
[[nodiscard]] auto HashStr32(std::string_view str) -> u32;

#endif  // SRC_LANCET_BASE_HASH_H_
