#include "lancet/hts/bgzf_ostream.h"

#include <cstdio>
#include <utility>

extern "C" {
#include "htslib/tbx.h"
}

namespace lancet::hts {

namespace detail {

auto BgzfStreambuf::Open(const std::filesystem::path &path, const char *mode) -> bool {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mFilePtr != nullptr) Close();

  mFileName = path;
  mFilePtr = bgzf_open(mFileName.c_str(), mode);
  return mFilePtr != nullptr;
}

void BgzfStreambuf::Close() {
  if (mFilePtr != nullptr) {
    bgzf_close(mFilePtr);
    mFilePtr = nullptr;
  }
}

auto BgzfStreambuf::uflow() -> int {
  if (mCurrPos != SENTINEL_BUFFER_POSITION) {
    const auto res = mCurrPos;
    mCurrPos = SENTINEL_BUFFER_POSITION;
    return res;
  }

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mFilePtr == nullptr) return EOF;

  mCurrPos = bgzf_getc(mFilePtr);
  switch (mCurrPos) {
    case -1:
    case -2:
      return EOF;
    default:
      return mCurrPos;
  }
}

auto BgzfStreambuf::underflow() -> int {
  // NOLINTBEGIN(readability-braces-around-statements)
  if (mFilePtr == nullptr) return EOF;
  if (mCurrPos != SENTINEL_BUFFER_POSITION) return mCurrPos;
  // NOLINTEND(readability-braces-around-statements)

  mCurrPos = bgzf_getc(mFilePtr);
  switch (mCurrPos) {
    case -1:
    case -2:
      return EOF;
    default:
      return mCurrPos;
  }
}

auto BgzfStreambuf::overflow(int dat) -> int {
  const auto cdat = static_cast<char>(dat);
  const auto num_bytes = bgzf_write(mFilePtr, &cdat, 1);
  return num_bytes < 0 ? static_cast<int>(num_bytes) : dat;
}

auto BgzfStreambuf::xsputn(const char *data, std::streamsize len) -> std::streamsize {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mFilePtr == nullptr) return 0;
  return bgzf_write(mFilePtr, data, static_cast<std::size_t>(len));
}

}  // namespace detail

auto BgzfOstream::Open(const std::filesystem::path &path, BgzfFormat ofmt) -> bool {
  mOutFmt = ofmt;
  auto result = mBgzfBuffer.Open(path, "w");
  rdbuf(&mBgzfBuffer);
  return result;
}

void BgzfOstream::Close() {
  mBgzfBuffer.Close();
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mOutFmt != BgzfFormat::UNSPECIFIED) BuildIndex();
}

void BgzfOstream::BuildIndex() {
  switch (mOutFmt) {
    case BgzfFormat::VCF:
      tbx_index_build(mBgzfBuffer.mFileName.c_str(), 0, &tbx_conf_vcf);
      break;
    case BgzfFormat::GFF:
      tbx_index_build(mBgzfBuffer.mFileName.c_str(), 0, &tbx_conf_gff);
      break;
    case BgzfFormat::BED:
      tbx_index_build(mBgzfBuffer.mFileName.c_str(), 0, &tbx_conf_bed);
      break;
    default:
      break;
  }
}

}  // namespace lancet::hts
