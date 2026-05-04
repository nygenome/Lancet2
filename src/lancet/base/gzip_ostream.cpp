#include "lancet/base/gzip_ostream.h"

#include "lancet/base/types.h"

#include "spdlog/fmt/bundled/format.h"
#include "zconf.h"
#include "zlib.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace lancet::base {

// 64 KiB output buffer balances syscall amortisation (one write() per
// drain) against working-set footprint (~64 KiB per active GzipOstream is
// negligible at any worker-pool size).
static constexpr usize DEFLATE_OUTPUT_BUFFER_BYTES = 65'536;

// windowBits = MAX_WBITS (15) selects the largest sliding window for best
// compression ratio; the +16 offset selects gzip wrapper format (writes
// the 10-byte gzip header + 8-byte trailer with CRC32 + ISIZE), as opposed
// to raw DEFLATE (offset 0) or zlib-format (offset 0 with no offset).
static constexpr int GZIP_WINDOW_BITS = MAX_WBITS + 16;

// Default deflate memory-level: 8 is zlib's documented default; covers
// ~256 KiB of internal compression state per stream — modest and matches
// what every other gzip-using tool uses.
static constexpr int DEFLATE_MEM_LEVEL = 8;

GzipOstream::GzipOstream(std::filesystem::path const& path, int const compression_level)
    : mOwnedFile(path, std::ios::binary | std::ios::trunc),
      mDeflateOutBuf(DEFLATE_OUTPUT_BUFFER_BYTES, u8{0}) {
  if (!mOwnedFile.is_open()) {
    throw std::runtime_error(
        fmt::format("GzipOstream: failed to open output file: {}", path.string()));
  }
  mSink = &mOwnedFile;
  InitDeflateStream(compression_level);
}

GzipOstream::GzipOstream(std::ostream& sink, int const compression_level)
    : mSink(&sink), mDeflateOutBuf(DEFLATE_OUTPUT_BUFFER_BYTES, u8{0}) {
  InitDeflateStream(compression_level);
}

void GzipOstream::InitDeflateStream(int const compression_level) {
  // mDeflateStream is value-initialised in the class body, which zeroes
  // every field (including zalloc/zfree/opaque, equivalent to Z_NULL).
  auto const init_result = deflateInit2(&mDeflateStream, compression_level, Z_DEFLATED,
                                        GZIP_WINDOW_BITS, DEFLATE_MEM_LEVEL, Z_DEFAULT_STRATEGY);
  if (init_result != Z_OK) {
    throw std::runtime_error(
        fmt::format("GzipOstream: deflateInit2 failed with code {}", init_result));
  }
}

GzipOstream::~GzipOstream() {
  // Best-effort close. If Close() throws here, swallow it: a destructor
  // throwing during stack unwinding terminates the program. Callers who
  // want to observe close errors should call Close() explicitly.
  try {
    Close();
    // Catch-all to avoid destructor throwing during stack unwinding
    // NOLINTNEXTLINE(bugprone-empty-catch)
  } catch (...) {}
}

void GzipOstream::Write(std::string_view input_bytes) {
  if (mIsClosed) {
    throw std::runtime_error("GzipOstream::Write called after Close()");
  }
  if (input_bytes.empty()) return;

  // Reinterpret the input view as the unsigned-byte pointer zlib expects.
  // The const_cast is required because zlib's `next_in` is `Bytef*` (non-
  // const) even though deflate() only reads the input bytes; safe in
  // practice and the standard idiom for zlib's C API.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast)
  mDeflateStream.next_in = const_cast<u8*>(reinterpret_cast<u8 const*>(input_bytes.data()));
  mDeflateStream.avail_in = static_cast<uInt>(input_bytes.size());

  DeflateLoopAndDrain(Z_NO_FLUSH);
}

void GzipOstream::Close() {
  if (mIsClosed) return;

  // Drive deflate to completion: Z_FINISH tells zlib to emit any buffered
  // input, the gzip trailer, and signal Z_STREAM_END. After this loop
  // returns, the stream is fully written.
  mDeflateStream.next_in = nullptr;
  mDeflateStream.avail_in = 0;
  DeflateLoopAndDrain(Z_FINISH);

  auto const end_result = deflateEnd(&mDeflateStream);
  // Z_DATA_ERROR from deflateEnd just means some discarded input remained
  // — a no-op for our streaming usage. Surface only true error codes.
  if (end_result != Z_OK && end_result != Z_DATA_ERROR) {
    throw std::runtime_error(
        fmt::format("GzipOstream::Close: deflateEnd failed with code {}", end_result));
  }

  // Path-ctor mode: we own the backing file — close it. Sink-ctor mode:
  // the caller owns the sink — only flush the gzip trailer through.
  if (mOwnedFile.is_open()) {
    mOwnedFile.close();
    if (mOwnedFile.fail() && !mOwnedFile.eof()) {
      throw std::runtime_error("GzipOstream::Close: ofstream close reported failure");
    }
  } else if (mSink != nullptr) {
    mSink->flush();
  }
  mIsClosed = true;
}

void GzipOstream::DeflateLoopAndDrain(int const flush_mode) {
  // Inner loop invariant: each pass refills `next_out`/`avail_out` from
  // the reusable buffer, calls deflate(), then drains the bytes deflate
  // produced into the ofstream. Continues while either (a) input remains,
  // or (b) flush_mode==Z_FINISH and we haven't yet seen Z_STREAM_END.
  while (true) {
    mDeflateStream.next_out = mDeflateOutBuf.data();
    mDeflateStream.avail_out = static_cast<uInt>(mDeflateOutBuf.size());

    auto const deflate_result = deflate(&mDeflateStream, flush_mode);
    if (deflate_result == Z_STREAM_ERROR) {
      throw std::runtime_error("GzipOstream: deflate returned Z_STREAM_ERROR");
    }

    auto const produced_byte_count = mDeflateOutBuf.size() - mDeflateStream.avail_out;
    if (produced_byte_count > 0) {
      DrainDeflateBuffer(produced_byte_count);
    }

    auto const all_input_consumed = mDeflateStream.avail_in == 0;
    auto const finishing_complete = flush_mode == Z_FINISH && deflate_result == Z_STREAM_END;
    auto const not_finishing_and_no_more_input = flush_mode != Z_FINISH && all_input_consumed;
    if (finishing_complete || not_finishing_and_no_more_input) return;

    // For Z_FINISH we keep looping until Z_STREAM_END even after input is
    // consumed; the trailer + final flushed bytes still need to come out.
    if (flush_mode == Z_FINISH &&
        deflate_result != Z_OK &&
        deflate_result != Z_STREAM_END &&
        deflate_result != Z_BUF_ERROR) {
      throw std::runtime_error(
          fmt::format("GzipOstream: deflate(Z_FINISH) failed with code {}", deflate_result));
    }
  }
}

void GzipOstream::DrainDeflateBuffer(usize const byte_count) {
  // ostream::write takes char const*; the deflate output buffer is u8 (zlib `Bytef`).
  // The two byte-pointer types are layout-compatible but distinct under strict aliasing.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* const buffer_as_chars = reinterpret_cast<char const*>(mDeflateOutBuf.data());
  mSink->write(buffer_as_chars, static_cast<std::streamsize>(byte_count));
  if (!mSink->good()) {
    throw std::runtime_error("GzipOstream: ostream::write failed mid-stream");
  }
}

}  // namespace lancet::base
