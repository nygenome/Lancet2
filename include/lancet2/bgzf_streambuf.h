#pragma once

#include <streambuf>
#include <string>

extern "C" {
#include "htslib/bgzf.h"
}

namespace lancet2::detail {
class BgzfStreambuf : public std::streambuf {
 public:
  std::string FileName;  // NOLINT

  BgzfStreambuf() = default;
  BgzfStreambuf(const BgzfStreambuf& other) = default;
  BgzfStreambuf(BgzfStreambuf&& other) = default;
  auto operator=(const BgzfStreambuf& other) -> BgzfStreambuf& = default;
  auto operator=(BgzfStreambuf&& other) -> BgzfStreambuf& = default;

  ~BgzfStreambuf() override { close(); }

  auto open(const std::string& path, const char* mode) -> bool;

  void close();

  auto uflow() -> int override;

  auto underflow() -> int override;

  auto overflow(int c = EOF) -> int override;  // NOLINT

  auto xsputn(const char* s, std::streamsize n) -> std::streamsize override;

 private:
  BGZF* fp = nullptr;
  int currPos = 0;
};
}  // namespace lancet2::detail
