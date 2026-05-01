#ifndef TESTS_BASE_TAR_GZ_TEST_HELPERS_H_
#define TESTS_BASE_TAR_GZ_TEST_HELPERS_H_

#include "lancet/base/types.h"

#include "zconf.h"
#include "zlib.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace lancet::tests {

// Read all bytes from a file into memory. Used by round-trip tests.
[[nodiscard]] inline auto ReadAllBytes(std::filesystem::path const& path) -> std::vector<u8> {
  std::ifstream in_stream(path, std::ios::binary);
  if (!in_stream.is_open()) {
    throw std::runtime_error("ReadAllBytes: failed to open " + path.string());
  }
  in_stream.seekg(0, std::ios::end);
  auto const file_byte_count = static_cast<usize>(in_stream.tellg());
  in_stream.seekg(0, std::ios::beg);

  std::vector<u8> file_bytes(file_byte_count, u8{0});
  if (file_byte_count > 0) {
    // istream::read takes char*; our buffer is u8 — layout-compatible but distinct types.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    in_stream.read(reinterpret_cast<char*>(file_bytes.data()),
                   static_cast<std::streamsize>(file_byte_count));
  }
  return file_bytes;
}

// Decompress a gzip-format byte buffer using zlib. Handles multi-member gzip
// streams (RFC 1952) by explicitly resetting the inflate state at each
// member boundary, since stock zlib does not auto-restart on Z_STREAM_END
// when more input remains. Returns the concatenation of all members'
// decompressed payloads.
[[nodiscard]] inline auto DecompressGzip(std::vector<u8> const& gz_bytes) -> std::vector<u8> {
  std::vector<u8> decompressed;
  if (gz_bytes.empty()) return decompressed;

  z_stream inflate_stream{};
  inflate_stream.zalloc = Z_NULL;
  inflate_stream.zfree = Z_NULL;
  inflate_stream.opaque = Z_NULL;
  // zlib's `next_in` is non-const `Bytef*` despite read-only semantics — standard idiom for the C
  // API.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
  inflate_stream.next_in = const_cast<u8*>(gz_bytes.data());
  inflate_stream.avail_in = static_cast<uInt>(gz_bytes.size());

  // windowBits = 15 + 32 enables auto-detection of gzip vs zlib format on
  // the *first* member; we re-arm this for each additional member below.
  static constexpr int GZIP_AUTO_DETECT_WINDOW_BITS = 15 + 32;
  if (inflateInit2(&inflate_stream, GZIP_AUTO_DETECT_WINDOW_BITS) != Z_OK) {
    throw std::runtime_error("DecompressGzip: inflateInit2 failed");
  }

  static constexpr usize INFLATE_OUT_CHUNK_BYTES = 64 * 1024;
  std::vector<u8> chunk_buf(INFLATE_OUT_CHUNK_BYTES, u8{0});
  while (true) {
    inflate_stream.next_out = chunk_buf.data();
    inflate_stream.avail_out = static_cast<uInt>(chunk_buf.size());

    auto const inflate_result = inflate(&inflate_stream, Z_NO_FLUSH);
    if (inflate_result != Z_OK && inflate_result != Z_STREAM_END) {
      inflateEnd(&inflate_stream);
      throw std::runtime_error("DecompressGzip: inflate failed with code " +
                               std::to_string(inflate_result));
    }

    auto const produced = chunk_buf.size() - inflate_stream.avail_out;
    decompressed.insert(decompressed.end(), chunk_buf.begin(),
                        chunk_buf.begin() + static_cast<i64>(produced));

    if (inflate_result == Z_STREAM_END) {
      if (inflate_stream.avail_in == 0) break;
      // More input remains: reset for the next gzip member.
      if (inflateReset(&inflate_stream) != Z_OK) {
        inflateEnd(&inflate_stream);
        throw std::runtime_error("DecompressGzip: inflateReset failed");
      }
      continue;
    }
    // inflate_result == Z_OK here. Keep looping while either the output
    // buffer is full (more decompressed bytes pending) or input remains.
    if (inflate_stream.avail_out == 0) continue;
    if (inflate_stream.avail_in == 0) break;
  }

  inflateEnd(&inflate_stream);
  return decompressed;
}

// One regular-file entry parsed from a USTAR archive.
struct TarEntry {
  std::string mEntryPath;
  std::string mContents;
};

// Parse a USTAR-formatted byte buffer into a sequence of regular-file entries.
// Stops at the first all-zero 512-byte header (TAR end-of-archive marker) or
// when the buffer is exhausted.
[[nodiscard]] inline auto ParseTarEntries(std::vector<u8> const& tar_bytes)
    -> std::vector<TarEntry> {
  static constexpr usize TAR_BLOCK_BYTES = 512;
  static constexpr usize TAR_NAME_FIELD_OFFSET = 0;
  static constexpr usize TAR_NAME_FIELD_BYTES = 100;
  static constexpr usize TAR_SIZE_FIELD_OFFSET = 124;
  static constexpr usize TAR_SIZE_FIELD_BYTES = 12;
  static constexpr usize TAR_MAGIC_FIELD_OFFSET = 257;
  static constexpr usize TAR_MAGIC_FIELD_BYTES = 6;
  static constexpr usize TAR_PREFIX_FIELD_OFFSET = 345;
  static constexpr usize TAR_PREFIX_FIELD_BYTES = 155;

  std::vector<TarEntry> parsed_entries;
  usize cursor = 0;
  while (cursor + TAR_BLOCK_BYTES <= tar_bytes.size()) {
    auto const* const header = tar_bytes.data() + cursor;

    // Detect end-of-archive: an all-zero header block.
    bool is_zero_block = true;
    for (usize byte_index = 0; byte_index < TAR_BLOCK_BYTES; ++byte_index) {
      if (header[byte_index] != 0) {
        is_zero_block = false;
        break;
      }
    }
    if (is_zero_block) break;

    // Validate USTAR magic so a non-USTAR header doesn't get misread.
    std::string_view const magic_field(
        reinterpret_cast<char const*>(header + TAR_MAGIC_FIELD_OFFSET), TAR_MAGIC_FIELD_BYTES);
    if (magic_field.substr(0, 5) != "ustar") {
      throw std::runtime_error("ParseTarEntries: bad USTAR magic at block offset " +
                               std::to_string(cursor));
    }

    // Decode name + prefix into a single archive path.
    auto const read_nul_string = [](u8 const* field_start, usize field_width) -> std::string {
      auto const* field_end = field_start;
      auto const* field_limit = field_start + field_width;
      while (field_end < field_limit && *field_end != 0) ++field_end;
      // tar header fields are NUL-padded ASCII; std::string ctor takes char const* over our u8
      // buffer.
      // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
      return std::string(reinterpret_cast<char const*>(field_start),
                         static_cast<usize>(field_end - field_start));
    };
    auto const name_part = read_nul_string(header + TAR_NAME_FIELD_OFFSET, TAR_NAME_FIELD_BYTES);
    auto const prefix_part =
        read_nul_string(header + TAR_PREFIX_FIELD_OFFSET, TAR_PREFIX_FIELD_BYTES);
    auto const entry_path = prefix_part.empty() ? name_part : (prefix_part + "/" + name_part);

    // Decode size from the octal ASCII size field.
    std::string size_str;
    for (usize byte_index = 0; byte_index < TAR_SIZE_FIELD_BYTES; ++byte_index) {
      char const size_char = static_cast<char>(header[TAR_SIZE_FIELD_OFFSET + byte_index]);
      if (size_char == '\0' || size_char == ' ') break;
      size_str.push_back(size_char);
    }
    usize const content_byte_count =
        size_str.empty() ? usize{0} : std::stoull(size_str, nullptr, 8);

    cursor += TAR_BLOCK_BYTES;
    if (cursor + content_byte_count > tar_bytes.size()) {
      throw std::runtime_error("ParseTarEntries: content runs past buffer end");
    }
    // string_view ctor takes char const*; tar payload buffer is u8 — layout-compatible bytes.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    auto const content_view = std::string_view(
        reinterpret_cast<char const*>(tar_bytes.data() + cursor), content_byte_count);
    parsed_entries.push_back(TarEntry{entry_path, std::string(content_view)});

    auto const padded_content_bytes =
        content_byte_count == 0
            ? usize{0}
            : ((content_byte_count + TAR_BLOCK_BYTES - 1) / TAR_BLOCK_BYTES) * TAR_BLOCK_BYTES;
    cursor += padded_content_bytes;
  }
  return parsed_entries;
}

// Verify that the trailing 1024 bytes of `tar_bytes` are all zero (the
// POSIX TAR end-of-archive marker = two consecutive 512-byte zero blocks).
[[nodiscard]] inline auto HasEndOfArchiveMarker(std::vector<u8> const& tar_bytes) -> bool {
  static constexpr usize TAR_END_MARKER_BYTES = 1024;
  if (tar_bytes.size() < TAR_END_MARKER_BYTES) return false;
  for (usize byte_index = tar_bytes.size() - TAR_END_MARKER_BYTES; byte_index < tar_bytes.size();
       ++byte_index) {
    if (tar_bytes[byte_index] != 0) return false;
  }
  return true;
}

// Allocate a fresh per-test scratch directory under temp_directory_path
// and ensure it's empty. Returns the directory path.
[[nodiscard]] inline auto MakeFreshScratchDir(std::string_view test_name) -> std::filesystem::path {
  auto const scratch_dir = std::filesystem::temp_directory_path() / std::string(test_name);
  std::filesystem::remove_all(scratch_dir);
  std::filesystem::create_directories(scratch_dir);
  return scratch_dir;
}

}  // namespace lancet::tests

#endif  // TESTS_BASE_TAR_GZ_TEST_HELPERS_H_
