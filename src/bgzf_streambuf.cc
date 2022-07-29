#include "lancet2/bgzf_streambuf.h"

#include <cstdio>
#include <utility>

namespace lancet2::detail {
auto BgzfStreambuf::open(const std::string &path, const char *mode) -> bool {
  if (fp != nullptr) close();

  FileName = path;
  fp = bgzf_open(FileName.c_str(), mode);
  return fp != nullptr;
}

void BgzfStreambuf::close() {
  if (fp != nullptr) {
    bgzf_close(fp);
    fp = nullptr;
  }
}

auto BgzfStreambuf::uflow() -> int {
  if (currPos != -999) {
    const auto res = currPos;
    currPos = -999;
    return res;
  }

  if (fp == nullptr) return EOF;

  currPos = bgzf_getc(fp);
  switch (currPos) {
    case -1:
    case -2:
      return EOF;
    default:
      return currPos;
  }
}

auto BgzfStreambuf::underflow() -> int {
  if (fp == nullptr) return EOF;
  if (currPos != -999) return currPos;

  currPos = bgzf_getc(fp);
  switch (currPos) {
    case -1:
    case -2:
      return EOF;
    default:
      return currPos;
  }
}

auto BgzfStreambuf::overflow(int c) -> int {  // NOLINT
  const auto z = static_cast<char>(c);
  const auto numBytes = bgzf_write(fp, &z, 1);
  return numBytes < 0 ? static_cast<int>(numBytes) : c;
}

auto BgzfStreambuf::xsputn(const char *s, std::streamsize n) -> std::streamsize {
  if (fp == nullptr) return 0;
  return bgzf_write(fp, s, static_cast<std::size_t>(n));
}
}  // namespace lancet2::detail
