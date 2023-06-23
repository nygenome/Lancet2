#ifndef SRC_LANCET_BASE_HASH_H_
#define SRC_LANCET_BASE_HASH_H_

#include <string_view>

#include "lancet/base/types.h"

[[nodiscard]] auto HashStr(std::string_view str) -> u64;

#endif  // SRC_LANCET_BASE_HASH_H_
