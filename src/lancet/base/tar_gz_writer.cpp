#include "lancet/base/tar_gz_writer.h"

#include "lancet/base/gzip_ostream.h"
#include "lancet/base/types.h"

#include "spdlog/fmt/bundled/format.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include <cstring>

namespace lancet::base {

// ── USTAR header field byte offsets and sizes ──────────────────────────
// See `man 5 tar` (USTAR section) for the canonical layout. Numbers here
// are intentionally explicit (no struct overlay) so the field encoding is
// readable and there's no platform-dependent padding.

static constexpr usize TAR_BLOCK_BYTES = 512;
static constexpr usize TAR_HEADER_BYTES = TAR_BLOCK_BYTES;
static constexpr usize TAR_NAME_FIELD_BYTES = 100;
static constexpr usize TAR_NAME_FIELD_OFFSET = 0;
static constexpr usize TAR_MODE_FIELD_BYTES = 8;
static constexpr usize TAR_MODE_FIELD_OFFSET = 100;
static constexpr usize TAR_UID_FIELD_BYTES = 8;
static constexpr usize TAR_UID_FIELD_OFFSET = 108;
static constexpr usize TAR_GID_FIELD_BYTES = 8;
static constexpr usize TAR_GID_FIELD_OFFSET = 116;
static constexpr usize TAR_SIZE_FIELD_BYTES = 12;
static constexpr usize TAR_SIZE_FIELD_OFFSET = 124;
static constexpr usize TAR_MTIME_FIELD_BYTES = 12;
static constexpr usize TAR_MTIME_FIELD_OFFSET = 136;
static constexpr usize TAR_CHKSUM_FIELD_BYTES = 8;
static constexpr usize TAR_CHKSUM_FIELD_OFFSET = 148;
static constexpr usize TAR_TYPEFLAG_FIELD_OFFSET = 156;
static constexpr usize TAR_MAGIC_FIELD_BYTES = 6;
static constexpr usize TAR_MAGIC_FIELD_OFFSET = 257;
static constexpr usize TAR_VERSION_FIELD_BYTES = 2;
static constexpr usize TAR_VERSION_FIELD_OFFSET = 263;
static constexpr usize TAR_PREFIX_FIELD_BYTES = 155;
static constexpr usize TAR_PREFIX_FIELD_OFFSET = 345;

// USTAR magic + version. The trailing NUL on "ustar" is part of the 6-byte
// MAGIC field; "00" is the literal version string (no NUL).
static constexpr std::string_view TAR_USTAR_MAGIC = "ustar";  // NUL-terminated in 6 bytes
static constexpr std::string_view TAR_USTAR_VERSION_DIGITS = "00";

// USTAR typeflag for a regular file. POSIX accepts both '0' and NUL ('\0')
// for back-compat with v7 tar; we use '0' to match GNU/BSD tar exactly.
static constexpr char TAR_TYPEFLAG_REGULAR_FILE = '0';

// Hardcoded mode/uid/gid/mtime: regular files at 0644, root:root, mtime
// of 0 (epoch). Using fixed values makes the archive deterministic and
// keeps Lancet's TAR writer free of any system-state dependence
// (tar -tvf will show "1970-01-01" — perfectly fine for debug artifacts).
static constexpr std::string_view TAR_DEFAULT_MODE_OCTAL = "000644";
static constexpr std::string_view TAR_DEFAULT_UID_OCTAL = "000000";
static constexpr std::string_view TAR_DEFAULT_GID_OCTAL = "000000";
static constexpr std::string_view TAR_DEFAULT_MTIME_OCTAL = "00000000000";

namespace {

// Reusable end-of-archive payload: two consecutive zero-filled 512-byte
// blocks. POSIX tar tooling stops at this marker (or warns "lone zero
// block at N" when only one is found). Constructed once at first use via
// the static `ZERO_BLOCKS` array below.
constexpr usize TAR_END_OF_ARCHIVE_BYTES = TAR_BLOCK_BYTES * 2;

[[nodiscard]] auto MakeZeroBlock() -> std::array<char, TAR_END_OF_ARCHIVE_BYTES> {
  std::array<char, TAR_END_OF_ARCHIVE_BYTES> zeros{};
  return zeros;
}

auto const& ZeroBlocks() {
  static auto const ZERO_BLOCKS = MakeZeroBlock();
  return ZERO_BLOCKS;
}

// Write `field_value` into `header_buf[offset .. offset + width)`, NUL-
// padding any trailing bytes. Throws if `field_value` exceeds `width`
// (callers ensure values fit).
void WriteFixedField(std::array<char, TAR_HEADER_BYTES>& header_buf, usize const offset,
                     usize const field_width, std::string_view const field_value) {
  if (field_value.size() > field_width) {
    throw std::runtime_error(fmt::format("TarGzWriter: USTAR field overflow at offset {}", offset));
  }
  std::memcpy(header_buf.data() + offset, field_value.data(), field_value.size());
  // Remaining bytes are already 0-initialised in the caller.
}

// Compute the USTAR header checksum: unsigned sum of all 512 bytes with
// the chksum field treated as 8 ASCII spaces. Result is rendered as
// "<6 octal digits>\0 " into the chksum field per the POSIX tar spec.
void WriteChecksumField(std::array<char, TAR_HEADER_BYTES>& header_buf) {
  // Fill the chksum slot with spaces for the sum calculation.
  std::memset(header_buf.data() + TAR_CHKSUM_FIELD_OFFSET, ' ', TAR_CHKSUM_FIELD_BYTES);

  u32 byte_sum = 0;
  for (auto const header_byte : header_buf) {
    byte_sum += static_cast<u8>(header_byte);
  }

  // POSIX format: 6 octal digits, NUL, space.
  auto const chksum_digits = fmt::format("{:06o}", byte_sum);
  std::memcpy(header_buf.data() + TAR_CHKSUM_FIELD_OFFSET, chksum_digits.data(),
              chksum_digits.size());
  header_buf[TAR_CHKSUM_FIELD_OFFSET + 6] = '\0';
  header_buf[TAR_CHKSUM_FIELD_OFFSET + 7] = ' ';
}

// Split `entry_path` between USTAR's prefix (155 chars) and name (100
// chars) fields at a `/` boundary. Throws if the path can't be split
// (>255 total, or no `/` in a position that respects both limits).
void SplitPathIntoNameAndPrefix(std::string_view const entry_path, std::string_view& out_prefix,
                                std::string_view& out_name) {
  if (entry_path.size() <= TAR_NAME_FIELD_BYTES) {
    out_prefix = {};
    out_name = entry_path;
    return;
  }

  if (entry_path.size() > TAR_NAME_FIELD_BYTES + TAR_PREFIX_FIELD_BYTES + 1) {
    throw std::runtime_error(
        fmt::format("TarGzWriter: entry path exceeds 255 chars (USTAR limit without LongLink): {}",
                    entry_path));
  }

  // Find the last '/' such that the name (everything after it) fits in
  // 100 chars AND the prefix (everything before it, exclusive) fits in
  // 155 chars. Iterate from the latest possible split toward the start.
  auto const max_split_index = std::min(entry_path.size() - 1, TAR_PREFIX_FIELD_BYTES);
  for (usize split_index = max_split_index; split_index > 0; --split_index) {
    if (entry_path[split_index] != '/') continue;

    auto const candidate_prefix_len = split_index;
    auto const candidate_name_len = entry_path.size() - split_index - 1;
    if (candidate_prefix_len <= TAR_PREFIX_FIELD_BYTES &&
        candidate_name_len <= TAR_NAME_FIELD_BYTES) {
      out_prefix = entry_path.substr(0, candidate_prefix_len);
      out_name = entry_path.substr(split_index + 1);
      return;
    }
  }

  static constexpr auto ERROR_MSG = "TarGzWriter: entry path > 100 chars cannot be split "
                                    "at a `/` boundary that respects USTAR's 100/155 split: {}";
  throw std::runtime_error(fmt::format(ERROR_MSG, entry_path));
}

// Build a 512-byte USTAR header for a regular-file entry of size
// `content_byte_count` at archive path `entry_path`.
[[nodiscard]] auto BuildUstarHeader(std::string_view const entry_path,
                                    usize const content_byte_count)
    -> std::array<char, TAR_HEADER_BYTES> {
  std::array<char, TAR_HEADER_BYTES> header_buf{};

  std::string_view path_prefix;
  std::string_view path_name;
  SplitPathIntoNameAndPrefix(entry_path, path_prefix, path_name);

  WriteFixedField(header_buf, TAR_NAME_FIELD_OFFSET, TAR_NAME_FIELD_BYTES, path_name);
  WriteFixedField(header_buf, TAR_PREFIX_FIELD_OFFSET, TAR_PREFIX_FIELD_BYTES, path_prefix);
  WriteFixedField(header_buf, TAR_MODE_FIELD_OFFSET, TAR_MODE_FIELD_BYTES, TAR_DEFAULT_MODE_OCTAL);
  WriteFixedField(header_buf, TAR_UID_FIELD_OFFSET, TAR_UID_FIELD_BYTES, TAR_DEFAULT_UID_OCTAL);
  WriteFixedField(header_buf, TAR_GID_FIELD_OFFSET, TAR_GID_FIELD_BYTES, TAR_DEFAULT_GID_OCTAL);
  WriteFixedField(header_buf, TAR_MTIME_FIELD_OFFSET, TAR_MTIME_FIELD_BYTES,
                  TAR_DEFAULT_MTIME_OCTAL);
  WriteFixedField(header_buf, TAR_MAGIC_FIELD_OFFSET, TAR_MAGIC_FIELD_BYTES, TAR_USTAR_MAGIC);
  WriteFixedField(header_buf, TAR_VERSION_FIELD_OFFSET, TAR_VERSION_FIELD_BYTES,
                  TAR_USTAR_VERSION_DIGITS);

  // Size: 11 octal digits + NUL, in a 12-byte field. Max representable
  // size is 8 GiB - 1; for any single Lancet entry this is irrelevant.
  auto const size_digits = fmt::format("{:011o}", content_byte_count);
  std::memcpy(header_buf.data() + TAR_SIZE_FIELD_OFFSET, size_digits.data(), size_digits.size());

  header_buf[TAR_TYPEFLAG_FIELD_OFFSET] = TAR_TYPEFLAG_REGULAR_FILE;

  // Checksum is computed last because it depends on all other fields.
  WriteChecksumField(header_buf);

  return header_buf;
}

}  // namespace

TarGzWriter::TarGzWriter(std::filesystem::path const& bundle_path,
                         EndOfArchive const end_marker_policy, int const compression_level)
    : mOutputSink(bundle_path, compression_level), mEndMarkerPolicy(end_marker_policy) {}

TarGzWriter::TarGzWriter(std::ostream& sink, EndOfArchive const end_marker_policy,
                         int const compression_level)
    : mOutputSink(sink, compression_level), mEndMarkerPolicy(end_marker_policy) {}

TarGzWriter::~TarGzWriter() {
  // Best-effort close. See GzipOstream::~GzipOstream for the same
  // rationale: throwing during stack unwinding terminates.
  try {
    Close();
    // intentional swallow: dtor cannot propagate exceptions
    // NOLINTNEXTLINE(bugprone-empty-catch)
  } catch (...) {}
}

void TarGzWriter::AddRegularFileEntry(std::string_view const entry_path,
                                      std::string_view const contents) {
  if (mIsClosed) {
    throw std::runtime_error("TarGzWriter::AddRegularFileEntry called after Close()");
  }

  auto const header_buf = BuildUstarHeader(entry_path, contents.size());
  mOutputSink.Write(std::string_view(header_buf.data(), header_buf.size()));

  if (!contents.empty()) {
    mOutputSink.Write(contents);
  }

  // Pad content + header up to a 512-byte boundary. The header is always
  // 512 bytes so the math is just on the contents length.
  auto const trailing_pad_bytes = contents.size() % TAR_BLOCK_BYTES == 0
                                      ? usize{0}
                                      : TAR_BLOCK_BYTES - (contents.size() % TAR_BLOCK_BYTES);
  if (trailing_pad_bytes > 0) {
    auto const& zeros = ZeroBlocks();
    mOutputSink.Write(std::string_view(zeros.data(), trailing_pad_bytes));
  }
}

void TarGzWriter::Close() {
  if (mIsClosed) return;

  if (mEndMarkerPolicy == EndOfArchive::EMIT) {
    auto const& zeros = ZeroBlocks();
    mOutputSink.Write(std::string_view(zeros.data(), TAR_END_OF_ARCHIVE_BYTES));
  }

  mOutputSink.Close();
  mIsClosed = true;
}

}  // namespace lancet::base
