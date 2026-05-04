#include "lancet/base/tar_gz_writer.h"

#include "lancet/base/types.h"

#include "catch_amalgamated.hpp"
#include "spdlog/fmt/bundled/format.h"
#include "tests/base/tar_gz_test_helpers.h"

#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lancet::base::tests {

using lancet::tests::DecompressGzip;
using lancet::tests::HasEndOfArchiveMarker;
using lancet::tests::MakeFreshScratchDir;
using lancet::tests::ParseTarEntries;
using lancet::tests::ReadAllBytes;

namespace {

// Helper: write one entry with EndOfArchive::EMIT, return the parsed
// archive listing for round-trip assertions.
[[nodiscard]] auto WriteSingleEntryAndParse(std::filesystem::path const& archive_path,
                                            std::string const& entry_path,
                                            std::string const& contents)
    -> std::vector<lancet::tests::TarEntry> {
  {
    TarGzWriter shard_writer(archive_path, TarGzWriter::EndOfArchive::EMIT);
    shard_writer.AddRegularFileEntry(entry_path, contents);
  }
  return ParseTarEntries(DecompressGzip(ReadAllBytes(archive_path)));
}

}  // namespace

TEST_CASE("TarGzWriter: name field exactly 100 chars (no prefix)", "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_name_100");
  auto const archive_path = scratch_dir / "name_100.tar.gz";

  auto const entry_path = std::string(100, 'a');
  auto const parsed_entries = WriteSingleEntryAndParse(archive_path, entry_path, "payload");
  REQUIRE(parsed_entries.size() == 1);
  CHECK(parsed_entries[0].mEntryPath == entry_path);
  CHECK(parsed_entries[0].mContents == "payload");

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter: path 101 chars splits into name + prefix", "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_split_101");
  auto const archive_path = scratch_dir / "split_101.tar.gz";

  // 50 chars + "/" + 50 chars = 101 chars. The split lands at the last `/`,
  // giving a 50-char prefix and a 50-char name (both fit comfortably).
  auto const prefix_part = std::string(50, 'p');
  auto const name_part = std::string(50, 'n');
  auto const entry_path = fmt::format("{}/{}", prefix_part, name_part);
  REQUIRE(entry_path.size() == 101);

  auto const parsed_entries = WriteSingleEntryAndParse(archive_path, entry_path, "split-101");
  REQUIRE(parsed_entries.size() == 1);
  CHECK(parsed_entries[0].mEntryPath == entry_path);
  CHECK(parsed_entries[0].mContents == "split-101");

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter: path near 255-char USTAR ceiling round-trips",
          "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_split_max");
  auto const archive_path = scratch_dir / "split_max.tar.gz";

  // 154-char prefix + "/" + 100-char name = 255 chars total. Exercises
  // both fields at their maximum width.
  auto const prefix_part = std::string(154, 'P');
  auto const name_part = std::string(100, 'N');
  auto const entry_path = fmt::format("{}/{}", prefix_part, name_part);
  REQUIRE(entry_path.size() == 255);

  auto const parsed_entries = WriteSingleEntryAndParse(archive_path, entry_path, "ceiling-payload");
  REQUIRE(parsed_entries.size() == 1);
  CHECK(parsed_entries[0].mEntryPath == entry_path);

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter: throws when entry path exceeds 255 chars", "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_too_long");
  auto const archive_path = scratch_dir / "too_long.tar.gz";

  // 156 + 100 = 256 chars exceeds the USTAR-without-LongLink ceiling.
  auto const overlong_entry_path =
      fmt::format("{}/{}", std::string(156, 'p'), std::string(99, 'n'));
  REQUIRE(overlong_entry_path.size() > 255);

  TarGzWriter shard_writer(archive_path, TarGzWriter::EndOfArchive::EMIT);
  CHECK_THROWS_AS(shard_writer.AddRegularFileEntry(overlong_entry_path, "x"), std::runtime_error);

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter: throws when path > 100 chars has no `/` for split",
          "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_no_slash");
  auto const archive_path = scratch_dir / "no_slash.tar.gz";

  // 150 chars, all letters, no `/` — there's no boundary to split on, so
  // the name field can't fit it on its own.
  auto const unsplittable_entry_path = std::string(150, 'q');

  TarGzWriter shard_writer(archive_path, TarGzWriter::EndOfArchive::EMIT);
  CHECK_THROWS_AS(shard_writer.AddRegularFileEntry(unsplittable_entry_path, "x"),
                  std::runtime_error);

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter: zero-byte entry round-trips", "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_empty_content");
  auto const archive_path = scratch_dir / "empty_content.tar.gz";

  auto const parsed_entries =
      WriteSingleEntryAndParse(archive_path, "dbg_graph/win_1/empty.dot", "");
  REQUIRE(parsed_entries.size() == 1);
  CHECK(parsed_entries[0].mEntryPath == "dbg_graph/win_1/empty.dot");
  CHECK(parsed_entries[0].mContents.empty());

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter: content size that's an exact 512-byte multiple skips padding",
          "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_aligned_size");
  auto const archive_path = scratch_dir / "aligned_size.tar.gz";

  // 1024 bytes of 'A' — exactly two TAR blocks, no trailing pad needed.
  auto const aligned_payload = std::string(1024, 'A');
  auto const parsed_entries =
      WriteSingleEntryAndParse(archive_path, "win/aligned.bin", aligned_payload);
  REQUIRE(parsed_entries.size() == 1);
  CHECK(parsed_entries[0].mContents == aligned_payload);

  // The decompressed bytes layout for one entry with 1024-byte content +
  // EOF marker = 512 (header) + 1024 (content, exact 2 blocks, no pad) +
  // 1024 (EOF) = 2560 bytes.
  auto const decompressed = DecompressGzip(ReadAllBytes(archive_path));
  static constexpr usize EXPECTED_BYTES = 512 + 1024 + 1024;
  CHECK(decompressed.size() == EXPECTED_BYTES);
  CHECK(HasEndOfArchiveMarker(decompressed));

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter: 1-byte content gets 511 bytes of trailing pad",
          "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_single_byte");
  auto const archive_path = scratch_dir / "single_byte.tar.gz";

  auto const parsed_entries = WriteSingleEntryAndParse(archive_path, "win/x.byte", "X");
  REQUIRE(parsed_entries.size() == 1);
  CHECK(parsed_entries[0].mContents == "X");

  // Layout: 512 (header) + 512 (1 content byte + 511 zero pad) + 1024 (EOF) = 2048.
  auto const decompressed = DecompressGzip(ReadAllBytes(archive_path));
  static constexpr usize EXPECTED_BYTES = 512 + 512 + 1024;
  CHECK(decompressed.size() == EXPECTED_BYTES);

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter: writes multiple entries in insertion order",
          "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_multi_entry");
  auto const archive_path = scratch_dir / "multi_entry.tar.gz";

  {
    TarGzWriter shard_writer(archive_path, TarGzWriter::EndOfArchive::EMIT);
    shard_writer.AddRegularFileEntry("dbg_graph/chr1_1_2/a.dot", "alpha");
    shard_writer.AddRegularFileEntry("dbg_graph/chr1_1_2/b.dot", "beta-content");
    shard_writer.AddRegularFileEntry("poa_graph/chr1_1_2/c.gfa", "gamma\nlines\n");
  }
  auto const parsed_entries = ParseTarEntries(DecompressGzip(ReadAllBytes(archive_path)));
  REQUIRE(parsed_entries.size() == 3);
  CHECK(parsed_entries[0].mEntryPath == "dbg_graph/chr1_1_2/a.dot");
  CHECK(parsed_entries[0].mContents == "alpha");
  CHECK(parsed_entries[1].mEntryPath == "dbg_graph/chr1_1_2/b.dot");
  CHECK(parsed_entries[1].mContents == "beta-content");
  CHECK(parsed_entries[2].mEntryPath == "poa_graph/chr1_1_2/c.gfa");
  CHECK(parsed_entries[2].mContents == "gamma\nlines\n");

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter: EndOfArchive::OMIT emits no trailing zero blocks",
          "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_omit_eof");
  auto const archive_path = scratch_dir / "omit_eof.tar.gz";

  {
    TarGzWriter shard_writer(archive_path, TarGzWriter::EndOfArchive::OMIT);
    shard_writer.AddRegularFileEntry("win/single.txt", "no-eof-marker");
  }
  auto const decompressed = DecompressGzip(ReadAllBytes(archive_path));

  // Layout: 512 (header) + 512 (13-byte content + 499-byte pad) = 1024.
  // No trailing 1024-byte zero EOF blocks, so total decompressed size is
  // exactly 1024 — not 2048.
  static constexpr usize EXPECTED_BYTES = 512 + 512;
  CHECK(decompressed.size() == EXPECTED_BYTES);
  // The last 1024 bytes here are the pad-rounded content, NOT zero blocks
  // (the trailing 499 bytes of pad are zero, but the leading 13 bytes are
  // the content "no-eof-marker"), so HasEndOfArchiveMarker reports false.
  CHECK_FALSE(HasEndOfArchiveMarker(decompressed));

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter: EndOfArchive::EMIT writes the 1024-byte zero trailer",
          "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_emit_eof");
  auto const archive_path = scratch_dir / "emit_eof.tar.gz";

  {
    TarGzWriter shard_writer(archive_path, TarGzWriter::EndOfArchive::EMIT);
    shard_writer.AddRegularFileEntry("win/single.txt", "with-eof-marker");
  }
  auto const decompressed = DecompressGzip(ReadAllBytes(archive_path));
  CHECK(HasEndOfArchiveMarker(decompressed));

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter::AddRegularFileEntry after Close throws", "[lancet][base][TarGzWriter]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_writer_post_close");
  auto const archive_path = scratch_dir / "post_close.tar.gz";

  TarGzWriter shard_writer(archive_path, TarGzWriter::EndOfArchive::EMIT);
  shard_writer.AddRegularFileEntry("win/before.txt", "before");
  shard_writer.Close();
  CHECK_THROWS_AS(shard_writer.AddRegularFileEntry("win/after.txt", "after"), std::runtime_error);

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzWriter sink ctor round-trips an entry through an in-memory stringstream",
          "[lancet][base][TarGzWriter]") {
  // The sink ctor delegates the std::ostream& to the internal GzipOstream,
  // letting tests assemble the tar.gz archive entirely in memory. Verifies
  // the same end-to-end shape (USTAR header + content + EOF marker) without
  // touching the filesystem.
  std::stringstream sink(std::ios::in | std::ios::out | std::ios::binary);
  {
    TarGzWriter shard_writer(sink, TarGzWriter::EndOfArchive::EMIT);
    shard_writer.AddRegularFileEntry("win/in_memory.txt", "no temp file needed");
  }

  auto const sink_str = sink.str();
  std::vector<u8> sink_bytes(sink_str.size(), u8{0});
  // string::data() returns char const*; our buffer is u8 — layout-compatible bytes.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* const src_bytes = reinterpret_cast<u8 const*>(sink_str.data());
  for (usize idx = 0; idx < sink_str.size(); ++idx) sink_bytes[idx] = src_bytes[idx];

  auto const decompressed = DecompressGzip(sink_bytes);
  auto const parsed_entries = ParseTarEntries(decompressed);
  REQUIRE(parsed_entries.size() == 1);
  CHECK(parsed_entries[0].mEntryPath == "win/in_memory.txt");
  CHECK(parsed_entries[0].mContents == "no temp file needed");
  CHECK(HasEndOfArchiveMarker(decompressed));
}

}  // namespace lancet::base::tests
