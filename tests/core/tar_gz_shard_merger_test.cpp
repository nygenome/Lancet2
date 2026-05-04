#include "lancet/core/tar_gz_shard_merger.h"

#include "lancet/base/tar_gz_writer.h"
#include "lancet/base/types.h"

#include "catch_amalgamated.hpp"
#include "spdlog/fmt/bundled/format.h"
#include "tests/base/tar_gz_test_helpers.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

using lancet::base::TarGzWriter;
using lancet::core::TarGzShardMerger;
using lancet::tests::DecompressGzip;
using lancet::tests::HasEndOfArchiveMarker;
using lancet::tests::MakeFreshScratchDir;
using lancet::tests::ParseTarEntries;
using lancet::tests::ReadAllBytes;

namespace {

// Helper: write a single-entry shard at `<shards_dir>/worker_<idx>.tar.gz`
// with EndOfArchive::OMIT (the layout the real workers produce).
void WriteWorkerShard(std::filesystem::path const& shards_dir, u32 const worker_index,
                      std::string const& entry_path, std::string const& contents) {
  auto const shard_path = shards_dir / fmt::format("worker_{}.tar.gz", worker_index);
  TarGzWriter shard_writer(shard_path, TarGzWriter::EndOfArchive::OMIT);
  shard_writer.AddRegularFileEntry(entry_path, contents);
}

}  // namespace

TEST_CASE("TarGzShardMerger: concatenates per-worker shards in worker-index order",
          "[lancet][core][TarGzShardMerger]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_shard_merger_basic");
  auto const shards_dir = scratch_dir / ".out.tar.gz.shards";
  auto const final_archive_path = scratch_dir / "out.tar.gz";
  std::filesystem::create_directories(shards_dir);

  // Note the deliberate out-of-order creation: the merger should still
  // process them in worker_<idx> ascending order regardless of fs
  // enumeration order.
  WriteWorkerShard(shards_dir, 2, "dbg_graph/win_C/c.dot", "from-worker-2");
  WriteWorkerShard(shards_dir, 0, "dbg_graph/win_A/a.dot", "from-worker-0");
  WriteWorkerShard(shards_dir, 1, "dbg_graph/win_B/b.dot", "from-worker-1");

  TarGzShardMerger merger(shards_dir, final_archive_path);
  merger.Merge();

  // Shards directory must be removed on success.
  CHECK_FALSE(std::filesystem::exists(shards_dir));

  // The final archive decompresses to <hdr><entry-0><pad><hdr><entry-1>
  // <pad><hdr><entry-2><pad><1024-zero EOF marker>.
  auto const decompressed = DecompressGzip(ReadAllBytes(final_archive_path));
  CHECK(HasEndOfArchiveMarker(decompressed));

  auto const parsed_entries = ParseTarEntries(decompressed);
  REQUIRE(parsed_entries.size() == 3);
  CHECK(parsed_entries[0].mEntryPath == "dbg_graph/win_A/a.dot");
  CHECK(parsed_entries[0].mContents == "from-worker-0");
  CHECK(parsed_entries[1].mEntryPath == "dbg_graph/win_B/b.dot");
  CHECK(parsed_entries[1].mContents == "from-worker-1");
  CHECK(parsed_entries[2].mEntryPath == "dbg_graph/win_C/c.dot");
  CHECK(parsed_entries[2].mContents == "from-worker-2");

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzShardMerger: parses worker indices as integers, not lexically",
          "[lancet][core][TarGzShardMerger]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_shard_merger_numeric_sort");
  auto const shards_dir = scratch_dir / ".out.tar.gz.shards";
  auto const final_archive_path = scratch_dir / "out.tar.gz";
  std::filesystem::create_directories(shards_dir);

  // 10, 2, 1 — under lexicographic ordering "10" < "2", but the merger
  // parses the integer index, so the order should be 1, 2, 10.
  WriteWorkerShard(shards_dir, 10, "win_X/a", "ten");
  WriteWorkerShard(shards_dir, 2, "win_Y/a", "two");
  WriteWorkerShard(shards_dir, 1, "win_Z/a", "one");

  TarGzShardMerger merger(shards_dir, final_archive_path);
  merger.Merge();

  auto const parsed_entries = ParseTarEntries(DecompressGzip(ReadAllBytes(final_archive_path)));
  REQUIRE(parsed_entries.size() == 3);
  CHECK(parsed_entries[0].mContents == "one");
  CHECK(parsed_entries[1].mContents == "two");
  CHECK(parsed_entries[2].mContents == "ten");

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzShardMerger: ignores files in the shards dir whose names don't match the pattern",
          "[lancet][core][TarGzShardMerger]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_shard_merger_unrelated_files");
  auto const shards_dir = scratch_dir / ".out.tar.gz.shards";
  auto const final_archive_path = scratch_dir / "out.tar.gz";
  std::filesystem::create_directories(shards_dir);

  WriteWorkerShard(shards_dir, 0, "win/a", "real-shard");

  // Drop unrelated files in the shards dir. The merger must skip them,
  // not include them in the merge stream.
  std::ofstream(shards_dir / "README.txt") << "human-readable note";
  std::ofstream(shards_dir / "worker_3.partial") << "incomplete shard";
  std::ofstream(shards_dir / ".hidden") << "hidden file";

  TarGzShardMerger merger(shards_dir, final_archive_path);
  merger.Merge();

  auto const parsed_entries = ParseTarEntries(DecompressGzip(ReadAllBytes(final_archive_path)));
  REQUIRE(parsed_entries.size() == 1);
  CHECK(parsed_entries[0].mContents == "real-shard");

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzShardMerger: produces an EOF-only archive when there are no shards",
          "[lancet][core][TarGzShardMerger]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_shard_merger_empty_shards");
  auto const shards_dir = scratch_dir / ".out.tar.gz.shards";
  auto const final_archive_path = scratch_dir / "out.tar.gz";
  std::filesystem::create_directories(shards_dir);

  TarGzShardMerger merger(shards_dir, final_archive_path);
  merger.Merge();

  auto const decompressed = DecompressGzip(ReadAllBytes(final_archive_path));
  // Just the 1024-byte EOF marker, nothing else.
  CHECK(decompressed.size() == 1024);
  CHECK(HasEndOfArchiveMarker(decompressed));
  CHECK(ParseTarEntries(decompressed).empty());

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzShardMerger: throws and preserves shards when final archive can't be opened",
          "[lancet][core][TarGzShardMerger]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_shard_merger_open_fail");
  auto const shards_dir = scratch_dir / ".out.tar.gz.shards";
  std::filesystem::create_directories(shards_dir);
  WriteWorkerShard(shards_dir, 0, "win/a", "preserve-me");

  // Point the final archive at a directory path that doesn't exist — the
  // ofstream open will fail and Merge() must throw without removing the
  // shards directory.
  auto const bogus_archive_path =
      scratch_dir / "no_such_subdir" / "no_such_other_subdir" / "out.tar.gz";
  TarGzShardMerger merger(shards_dir, bogus_archive_path);
  CHECK_THROWS_AS(merger.Merge(), std::runtime_error);
  CHECK(std::filesystem::exists(shards_dir));
  CHECK(std::filesystem::exists(shards_dir / "worker_0.tar.gz"));

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("TarGzShardMerger: multi-entry shards round-trip end-to-end",
          "[lancet][core][TarGzShardMerger]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_tar_gz_shard_merger_multi_entry_shards");
  auto const shards_dir = scratch_dir / ".out.tar.gz.shards";
  auto const final_archive_path = scratch_dir / "out.tar.gz";
  std::filesystem::create_directories(shards_dir);

  // Worker 0 writes two entries; worker 1 writes one. Both shards must
  // omit the EOF marker per the OMIT contract; the merger appends a
  // single EOF marker at the very end.
  {
    TarGzWriter shard_writer_0(shards_dir / "worker_0.tar.gz", TarGzWriter::EndOfArchive::OMIT);
    shard_writer_0.AddRegularFileEntry("dbg_graph/w0a/x.dot", "w0-entry-x");
    shard_writer_0.AddRegularFileEntry("poa_graph/w0a/y.gfa", "w0-entry-y");
  }
  {
    TarGzWriter shard_writer_1(shards_dir / "worker_1.tar.gz", TarGzWriter::EndOfArchive::OMIT);
    shard_writer_1.AddRegularFileEntry("dbg_graph/w1a/z.dot", "w1-entry-z");
  }

  TarGzShardMerger merger(shards_dir, final_archive_path);
  merger.Merge();

  auto const parsed_entries = ParseTarEntries(DecompressGzip(ReadAllBytes(final_archive_path)));
  REQUIRE(parsed_entries.size() == 3);
  CHECK(parsed_entries[0].mContents == "w0-entry-x");
  CHECK(parsed_entries[1].mContents == "w0-entry-y");
  CHECK(parsed_entries[2].mContents == "w1-entry-z");

  std::filesystem::remove_all(scratch_dir);
}
