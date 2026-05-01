#ifndef SRC_LANCET_HTS_BGZF_OSTREAM_H_
#define SRC_LANCET_HTS_BGZF_OSTREAM_H_

#include "lancet/base/types.h"

extern "C" {
#include "htslib/bgzf.h"
}

#include <filesystem>
#include <ios>
#include <ostream>
#include <stdio.h>

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

  ~BgzfStreambuf() override {
    // Destructors are implicitly noexcept; Close() may allocate when formatting
    // its error-path log message. Swallow any propagated exception so unwind
    // doesn't trigger std::terminate. Explicit Close() callers still see throws.
    try {
      Close();
      // intentional swallow: dtor cannot propagate exceptions
      // NOLINTNEXTLINE(bugprone-empty-catch)
    } catch (...) {}
  }

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
  // ── 8B Align ────────────────────────────────────────────────────────────
  BGZF* mFilePtr = nullptr;
  // ── 4B Align ────────────────────────────────────────────────────────────
  int mCurrPos = 0;
};

}  // namespace detail

enum class BgzfFormat : u8 { UNSPECIFIED, GFF, BED, VCF };

/// Thread safety: NOT thread-safe. Each writer thread must use its own
/// BgzfOstream instance. Concurrent writes to the same underlying BGZF
/// handle corrupt the block structure.
// NOLINTNEXTLINE(misc-multiple-inheritance)
class BgzfOstream : public std::ostream {
 public:
  BgzfOstream() : std::ostream(nullptr) {}
  ~BgzfOstream() override {
    // See BgzfStreambuf::~BgzfStreambuf for rationale.
    try {
      Close();
      // intentional swallow: dtor cannot propagate exceptions
      // NOLINTNEXTLINE(bugprone-empty-catch)
    } catch (...) {}
  }

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
  // ── 8B Align ────────────────────────────────────────────────────────────
  detail::BgzfStreambuf mBgzfBuffer;
  // ── 1B Align ────────────────────────────────────────────────────────────
  BgzfFormat mOutFmt = BgzfFormat::UNSPECIFIED;

  void BuildIndex();
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_BGZF_OSTREAM_H_
