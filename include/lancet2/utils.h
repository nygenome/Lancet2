#pragma once

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "absl/status/status.h"
#include "absl/types/span.h"
#include "lancet2/sized_ints.h"

namespace lancet2::utils {
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-const-variable"
#endif
static constexpr u64 PRIME_0 = 18446744073709551557LLU;
static constexpr u64 PRIME_1 = 14480561146010017169LLU;
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

[[nodiscard]] auto HammingDistWithin(std::string_view s1, std::string_view s2, usize max) -> bool;

[[nodiscard]] auto HasRepeatKmer(std::string_view seq, usize k) -> bool;

[[nodiscard]] auto HasAlmostRepeatKmer(std::string_view seq, usize k, u32 max_mismatches) -> bool;

[[nodiscard]] auto RevComp(const char& b) -> char;

[[nodiscard]] auto RevComp(std::string_view sv) -> std::string;

[[nodiscard]] auto RevStr(std::string_view sv) -> std::string;

[[nodiscard]] auto GetHash(std::string_view sv) -> u64;

void PushSeq(std::string_view sequence, std::string* result);

void PushRevCompSeq(std::string_view sequence, std::string* result);

void PushSeq(std::string* result, std::string::iterator position, std::string_view sequence, usize start_offset,
             usize end_offset);

void PushRevCompSeq(std::string* result, std::string::iterator position, std::string_view sequence,
                    usize rev_start_offset, usize rev_end_offset);

}  // namespace lancet2::utils
