#ifndef SRC_LANCET_HTS_BGZF_OSTREAM_H_
#define SRC_LANCET_HTS_BGZF_OSTREAM_H_

#include <filesystem>
#include <ios>
#include <ostream>
#include <streambuf>

extern "C" {
#include "htslib/bgzf.h"
}

#include "lancet/base/types.h"

namespace lancet::hts {

namespace detail {

class BgzfStreambuf : public std::streambuf {
 public:
  std::filesystem::path mFileName{};

  BgzfStreambuf() = default;
  BgzfStreambuf(BgzfStreambuf const& other) = default;
  BgzfStreambuf(BgzfStreambuf&& other) = default;
  auto operator=(BgzfStreambuf const& other) -> BgzfStreambuf& = default;
  auto operator=(BgzfStreambuf&& other) -> BgzfStreambuf& = default;

  ~BgzfStreambuf() override { Close(); }

  auto Open(std::filesystem::path const& path, char const* mode) -> bool;
  void Close();

  // NOLINTBEGIN(misc-override-with-different-visibility)
  auto uflow() -> int override;
  auto underflow() -> int override;
  auto overflow(int dat = EOF) -> int override;
  auto xsputn(char const* data, std::streamsize len) -> std::streamsize override;
  auto sync() -> int override;
  // NOLINTEND(misc-override-with-different-visibility)

 private:
  static constexpr int SENTINEL_BUFFER_POSITION = -999;
  BGZF* mFilePtr = nullptr;
  int mCurrPos = 0;
};

}  // namespace detail

enum class BgzfFormat : u8 { UNSPECIFIED, GFF, BED, VCF };

// NOLINTNEXTLINE(misc-multiple-inheritance)
class BgzfOstream : public std::ostream {
 public:
  BgzfOstream() : std::ostream(nullptr) {}
  ~BgzfOstream() override { Close(); }

  BgzfOstream(BgzfOstream const&) = delete;
  BgzfOstream(BgzfOstream&&) = delete;
  auto operator=(BgzfOstream const&) -> BgzfOstream& = delete;
  auto operator=(BgzfOstream&&) -> BgzfOstream& = delete;

  auto Open(std::filesystem::path const& path, BgzfFormat ofmt) -> bool;
  auto Open(std::filesystem::path const& path) -> bool {
    return Open(path, BgzfFormat::UNSPECIFIED);
  }
  void Close();

 private:
  detail::BgzfStreambuf mBgzfBuffer;
  BgzfFormat mOutFmt = BgzfFormat::UNSPECIFIED;

  void BuildIndex();
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_BGZF_OSTREAM_H_
