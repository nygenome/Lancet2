#ifndef SRC_LANCET_HTS_BGZF_OSTREAM_H_
#define SRC_LANCET_HTS_BGZF_OSTREAM_H_

#include <filesystem>
#include <ostream>
#include <streambuf>
#include <string>

extern "C" {
#include "htslib/bgzf.h"
}

namespace lancet::hts {

namespace detail {

class BgzfStreambuf : public std::streambuf {
 public:
  // NOLINTBEGIN(misc-non-private-member-variables-in-classes,cppcoreguidelines-non-private-member-variables-in-classes)
  std::filesystem::path mFileName;
  // NOLINTEND(misc-non-private-member-variables-in-classes,cppcoreguidelines-non-private-member-variables-in-classes)

  BgzfStreambuf() = default;
  BgzfStreambuf(const BgzfStreambuf& other) = default;
  BgzfStreambuf(BgzfStreambuf&& other) = default;
  auto operator=(const BgzfStreambuf& other) -> BgzfStreambuf& = default;
  auto operator=(BgzfStreambuf&& other) -> BgzfStreambuf& = default;

  ~BgzfStreambuf() override { Close(); }

  auto Open(const std::filesystem::path& path, const char* mode) -> bool;
  void Close();

  auto uflow() -> int override;

  auto underflow() -> int override;

  auto overflow(int dat = EOF) -> int override;  // NOLINT

  auto xsputn(const char* data, std::streamsize len) -> std::streamsize override;

 private:
  static constexpr int SENTINEL_BUFFER_POSITION = -999;
  BGZF* mFilePtr = nullptr;
  int mCurrPos = 0;
};

}  // namespace detail

enum class BgzfFormat { UNSPECIFIED, GFF, BED, VCF };

class BgzfOstream : public std::ostream {
 public:
  BgzfOstream() : std::ostream(nullptr) {}
  ~BgzfOstream() override { Close(); }

  BgzfOstream(const BgzfOstream&) = delete;
  BgzfOstream(BgzfOstream&&) = delete;
  auto operator=(const BgzfOstream&) -> BgzfOstream& = delete;
  auto operator=(BgzfOstream&&) -> BgzfOstream& = delete;

  auto Open(const std::filesystem::path& path, BgzfFormat ofmt) -> bool;
  auto Open(const std::filesystem::path& path) -> bool { return Open(path, BgzfFormat::UNSPECIFIED); }
  void Close();

 private:
  detail::BgzfStreambuf mBgzfBuffer;
  BgzfFormat mOutFmt = BgzfFormat::UNSPECIFIED;

  void BuildIndex();
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_BGZF_OSTREAM_H_
