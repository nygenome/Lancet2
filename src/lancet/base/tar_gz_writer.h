#ifndef SRC_LANCET_BASE_TAR_GZ_WRITER_H_
#define SRC_LANCET_BASE_TAR_GZ_WRITER_H_

#include "lancet/base/gzip_ostream.h"
#include "lancet/base/types.h"

#include <filesystem>
#include <ostream>
#include <string_view>

namespace lancet::base {

// ============================================================================
// TarGzWriter — minimal POSIX USTAR-format writer that streams entries
// through a GzipOstream sink, producing a gzip-compressed TAR archive.
//
// Single-thread per instance. Each `AddRegularFileEntry` call writes one
// 512-byte USTAR header + file contents padded up to a 512-byte boundary
// into the open gzip stream. `Close()` (or the destructor) finalises the
// archive: it optionally writes the two trailing 512-byte zero blocks
// (TAR end-of-archive marker) and then flushes the gzip trailer.
//
// ── Why the EndOfArchive policy ─────────────────────────────────────────
// When workers each write a per-shard archive that is later cat-merged
// into a single final archive, individual shards must NOT terminate
// themselves with the end-of-archive zero blocks. The merge appends a
// single end-of-archive marker once. Use `EndOfArchive::OMIT` for shards
// and `EndOfArchive::EMIT` for standalone archives.
//
// ── Path length limits (USTAR) ──────────────────────────────────────────
// USTAR splits long paths between a 100-char `name` field and a 155-char
// `prefix` field, joined by a `/`. Total ≤ 255 chars works without
// extensions. `AddRegularFileEntry` throws if a path can't be represented
// (e.g. > 255 chars, or no `/` boundary that fits the 100/155 split).
// For Lancet's actual outputs this never trips: longest observed entry
// path is ~86 chars, GRCh38 worst-case decoy contigs ~120-130 chars.
// ============================================================================
class TarGzWriter {
 public:
  enum class EndOfArchive : u8 { OMIT, EMIT };

  /// Open a gzipped TAR file at `bundle_path` for writing. Truncates if
  /// exists. Throws on file open or zlib init failure.
  explicit TarGzWriter(std::filesystem::path const& bundle_path, EndOfArchive end_marker_policy,
                       int compression_level = GzipOstream::DEFAULT_COMPRESSION_LEVEL);

  /// Sink ctor: write the tar.gz archive to a caller-owned `std::ostream`.
  /// Delegates the sink to the internal `GzipOstream`. Useful for tests
  /// that want to inspect bytes without touching the filesystem. The sink
  /// must outlive the TarGzWriter.
  TarGzWriter(std::ostream& sink, EndOfArchive end_marker_policy,
              int compression_level = GzipOstream::DEFAULT_COMPRESSION_LEVEL);

  TarGzWriter(TarGzWriter const&) = delete;
  auto operator=(TarGzWriter const&) -> TarGzWriter& = delete;
  TarGzWriter(TarGzWriter&&) noexcept = delete;
  auto operator=(TarGzWriter&&) noexcept -> TarGzWriter& = delete;

  /// Calls Close() if not already.
  ~TarGzWriter();

  /// Append one regular-file entry. `entry_path` is the path stored in
  /// the archive; `contents` is the file body. Throws on path-too-long
  /// (>255 chars or unsplittable) or I/O error.
  void AddRegularFileEntry(std::string_view entry_path, std::string_view contents);

  /// Finalise the archive: write the end-of-archive marker (if EMIT was
  /// configured), flush the gzip stream, close the file. Idempotent.
  void Close();

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  GzipOstream mOutputSink;

  // ── 1B Align ────────────────────────────────────────────────────────────
  EndOfArchive mEndMarkerPolicy;
  bool mIsClosed = false;
};

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_TAR_GZ_WRITER_H_
