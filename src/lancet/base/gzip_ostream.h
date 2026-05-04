#ifndef SRC_LANCET_BASE_GZIP_OSTREAM_H_
#define SRC_LANCET_BASE_GZIP_OSTREAM_H_

#include "lancet/base/types.h"

#include "zlib.h"

#include <filesystem>
#include <fstream>
#include <string_view>
#include <vector>

namespace lancet::base {

// ============================================================================
// GzipOstream — streaming gzip writer wrapping zlib-ng's deflate() API.
//
// Single-thread per instance. Writes input bytes incrementally through a
// `z_stream` configured to emit the gzip wire format (header + DEFLATE
// payload + trailer). The underlying file descriptor is kept open for the
// lifetime of the GzipOstream; data flows: caller -> deflate() -> internal
// output buffer -> the std::ofstream backing file.
//
// Use `Write()` repeatedly to feed input in arbitrary-sized chunks; call
// `Close()` (or let the destructor handle it) to finalise the stream
// (deflate(Z_FINISH) writes the gzip trailer + CRC32) and close the fd.
//
// `Close()` is idempotent. Throws on zlib errors (Z_STREAM_ERROR /
// Z_MEM_ERROR / Z_BUF_ERROR misuse) or filesystem I/O failures so callers
// don't silently lose data.
// ============================================================================
class GzipOstream {
 public:
  // Default fast-compression level. Empirically gzip-1 hits ~250-400 MB/s
  // single-thread on modern CPUs while still cutting text DOT/GFA/FASTA
  // output by ~2-4x. Higher levels add CPU time without proportional ratio
  // gains for our workload.
  static constexpr int DEFAULT_COMPRESSION_LEVEL = 1;

  /// Open `path` for writing (truncates if exists). Initialises a zlib-ng
  /// `z_stream` with `windowBits = 15 + 16` so deflate emits gzip-format
  /// output (header + trailer auto-managed). Throws on file open failure
  /// or zlib init error. Internally opens a `std::ofstream` and forwards
  /// to the sink ctor below.
  explicit GzipOstream(std::filesystem::path const& path,
                       int compression_level = DEFAULT_COMPRESSION_LEVEL);

  /// Sink ctor: write the gzip stream to a caller-owned `std::ostream`.
  /// Useful for tests (write into a `std::stringstream` and inspect the
  /// bytes without touching the filesystem) and for any caller that already
  /// has an open ostream. The sink must outlive the GzipOstream. Close()
  /// flushes the gzip trailer through the sink but does NOT close the sink
  /// itself — that lifetime stays with the caller.
  explicit GzipOstream(std::ostream& sink, int compression_level = DEFAULT_COMPRESSION_LEVEL);

  GzipOstream(GzipOstream const&) = delete;
  auto operator=(GzipOstream const&) -> GzipOstream& = delete;
  GzipOstream(GzipOstream&&) noexcept = delete;
  auto operator=(GzipOstream&&) noexcept -> GzipOstream& = delete;

  /// Calls Close() if not already.
  ~GzipOstream();

  /// Feed `input_bytes` into the deflate stream. Compressed output is
  /// drained from the internal buffer into the file as the buffer fills.
  /// Throws on Z_STREAM_ERROR / I/O failure.
  void Write(std::string_view input_bytes);

  /// Flush remaining input via deflate(Z_FINISH), write the gzip trailer,
  /// close the file. Idempotent — calling Close() twice is a no-op.
  void Close();

 private:
  /// Initialise the zlib `z_stream` for gzip output at `compression_level`.
  /// Called from both ctors after the sink is wired up.
  void InitDeflateStream(int compression_level);

  /// Drive the deflate inner loop with `flush_mode` (typically Z_NO_FLUSH
  /// during writes, Z_FINISH on close). Drains the deflate output buffer
  /// to the backing sink as it fills. Returns when zlib has consumed all
  /// pending input AND (for Z_FINISH) has emitted Z_STREAM_END.
  void DeflateLoopAndDrain(int flush_mode);

  /// Helper: write `byte_count` bytes from `mDeflateOutBuf` to the
  /// underlying sink (via `mSink`).
  void DrainDeflateBuffer(usize byte_count);

  // ── 8B Align ────────────────────────────────────────────────────────────
  // Path-ctor mode: `mOwnedFile` is the open backing ofstream and `mSink`
  // points at it. Sink-ctor mode: `mOwnedFile` stays unopened; `mSink`
  // points at the caller's ostream. Close() inspects `mOwnedFile.is_open()`
  // to decide whether to close (path ctor) or just flush (sink ctor).
  std::ofstream mOwnedFile;
  std::ostream* mSink = nullptr;
  z_stream mDeflateStream{};
  std::vector<u8> mDeflateOutBuf;  // 64 KiB reusable scratch

  // ── 1B Align ────────────────────────────────────────────────────────────
  bool mIsClosed = false;
};

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_GZIP_OSTREAM_H_
