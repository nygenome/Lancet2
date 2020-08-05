#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "absl/status/status.h"
#include "absl/types/span.h"

namespace lancet::utils {
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-const-variable"
#endif
static constexpr std::uint64_t PRIME_0 = 18446744073709551557LLU;
static constexpr std::uint64_t PRIME_1 = 14480561146010017169LLU;
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

[[nodiscard]] auto MakeDir(const std::filesystem::path& dirname) -> absl::Status;

[[nodiscard]] auto HammingDistWithin(std::string_view s1, std::string_view s2, std::size_t max) -> bool;

[[nodiscard]] auto HasRepeatKmer(std::string_view seq, std::size_t k) -> bool;

[[nodiscard]] auto HasAlmostRepeatKmer(std::string_view seq, std::size_t k, std::uint32_t max_mismatches) -> bool;

[[nodiscard]] auto RevComp(const char& b) -> char;

[[nodiscard]] auto RevComp(std::string_view sv) -> std::string;

[[nodiscard]] auto RevStr(std::string_view sv) -> std::string;

template <typename T>
[[nodiscard]] auto MakeVector(absl::Span<const T> s) -> std::vector<T> {
  std::vector<T> result;
  result.reserve(s.size());
  result.insert(result.begin(), s.cbegin(), s.cend());
  return result;
}

template <typename T>
[[nodiscard]] auto ReverseVector(absl::Span<const T> s) -> std::vector<T> {
  std::vector<T> result;
  result.reserve(s.size());
  result.insert(result.begin(), s.crbegin(), s.crend());
  return result;
}

void PushSeq(std::string_view sequence, std::string* result);

void PushRevCompSeq(std::string_view sequence, std::string* result);

void PushSeq(std::string* result, std::string::iterator position, std::string_view sequence, std::size_t start_offset,
             std::size_t end_offset);

void PushRevCompSeq(std::string* result, std::string::iterator position, std::string_view sequence,
                    std::size_t rev_start_offset, std::size_t rev_end_offset);

}  // namespace lancet::utils
