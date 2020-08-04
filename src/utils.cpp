#include "lancet/utils.h"

#include <cerrno>
#include <cstring>

#include "absl/container/flat_hash_set.h"
#include "sys/stat.h"

namespace lancet::utils {
auto MakeDir(const std::filesystem::path& dirname) -> absl::Status {
  const auto status = mkdir(dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);  // NOLINT
  if (status == 0) return absl::OkStatus();
  return absl::InternalError(std::strerror(errno));
}

auto HammingDistWithin(std::string_view s1, std::string_view s2, std::size_t max) -> bool {
  if (s1.length() != s2.length()) return false;

  std::size_t differences = 0;
  for (std::size_t idx = 0; idx < s1.length(); idx++) {
    if (s1[idx] != s2[idx]) ++differences;
    if (differences > max) return false;
  }

  return true;
}

auto HasRepeatKmer(std::string_view seq, std::size_t k) -> bool {
  absl::flat_hash_set<std::string_view> mers;
  const auto endOffset = seq.length() - k + 1;
  mers.reserve(endOffset);
  for (std::size_t offset = 0; offset < endOffset; ++offset) {
    const auto itr = mers.emplace(absl::ClippedSubstr(seq, offset, k));
    if (!itr.second) return true;
  }
  return false;
}

auto HasAlmostRepeatKmer(std::string_view seq, std::size_t k, std::uint32_t max_mismatches) -> bool {
  const auto endOffset = seq.length() - k + 1;
  for (std::size_t offset1 = 0; offset1 < endOffset; ++offset1) {
    for (std::size_t offset2 = offset1 + k; offset2 < endOffset; ++offset2) {
      std::uint32_t numSeenMismatches = 0;

      for (std::size_t idx = 0; idx < k; idx++) {
        if (seq[offset1 + idx] != seq[offset2 + idx]) numSeenMismatches++;
        if (numSeenMismatches > max_mismatches) break;
      }

      if (numSeenMismatches <= max_mismatches) return true;
    }
  }
  return false;
}

auto RevComp(const char& b) -> char {
  switch (b) {
    case 'A':
    case 'a':
      return 'T';

    case 'C':
    case 'c':
      return 'G';

    case 'G':
    case 'g':
      return 'C';

    case 'T':
    case 't':
      return 'A';

    default:
      return 'N';
  }
}

auto RevComp(std::string_view sv) -> std::string {
  std::string result;
  result.reserve(sv.length());
  std::for_each(sv.crbegin(), sv.crend(), [&result](const char& b) { result.push_back(RevComp(b)); });
  return result;
}

auto RevStr(std::string_view sv) -> std::string {
  std::string result;
  result.reserve(sv.length());
  result.insert(result.begin(), sv.crbegin(), sv.crend());
  return result;
}

void PushSeq(std::string_view sequence, std::string* result) {
  result->insert(result->end(), sequence.cbegin(), sequence.cend());
}

void PushRevCompSeq(std::string_view sequence, std::string* result) {
  std::for_each(sequence.crbegin(), sequence.crend(),
                [&result](const char& b) { result->push_back(utils::RevComp(b)); });
}

void PushSeq(std::string* result, std::string::iterator position, std::string_view sequence, std::size_t start_offset,
             std::size_t end_offset) {
  result->insert(position, sequence.cbegin() + start_offset, sequence.cbegin() + end_offset);
}

void PushRevCompSeq(std::string* result, std::string::iterator position, std::string_view sequence,
                    std::size_t rev_start_offset, std::size_t rev_end_offset) {
  std::string tmpChunk(sequence.crbegin() + rev_start_offset, sequence.crbegin() + rev_end_offset);
  std::transform(tmpChunk.begin(), tmpChunk.end(), tmpChunk.begin(), [](const char& b) { return utils::RevComp(b); });
  result->insert(position, tmpChunk.cbegin(), tmpChunk.cend());
}
}  // namespace lancet::utils
