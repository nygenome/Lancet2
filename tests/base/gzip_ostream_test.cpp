#include "lancet/base/gzip_ostream.h"

#include "lancet/base/types.h"

#include "catch_amalgamated.hpp"
#include "tests/base/tar_gz_test_helpers.h"

#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using lancet::base::GzipOstream;
using lancet::tests::DecompressGzip;
using lancet::tests::MakeFreshScratchDir;
using lancet::tests::ReadAllBytes;

namespace {

// Helper: write `payload` through a fresh GzipOstream, return the
// decompressed bytes for round-trip comparison.
[[nodiscard]] auto RoundTripGzip(std::filesystem::path const& gz_path, std::string const& payload)
    -> std::string {
  {
    GzipOstream gz_out(gz_path);
    gz_out.Write(payload);
  }  // ~GzipOstream calls Close() — gzip trailer + file close happens here.
  auto const decompressed = DecompressGzip(ReadAllBytes(gz_path));
  // std::string ctor takes char const*; decompressed buffer is u8 (zlib byte type).
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* const decompressed_chars = reinterpret_cast<char const*>(decompressed.data());
  return std::string{decompressed_chars, decompressed.size()};
}

}  // namespace

TEST_CASE("GzipOstream: round-trips a known short string", "[lancet][base][GzipOstream]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_gzip_ostream_short");
  auto const gz_path = scratch_dir / "hello.gz";

  auto const payload = std::string("Hello, gzip world! Lancet2 streaming gzip writer test.\n");
  CHECK(RoundTripGzip(gz_path, payload) == payload);

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("GzipOstream: handles empty stream — gzip header + trailer only",
          "[lancet][base][GzipOstream]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_gzip_ostream_empty");
  auto const gz_path = scratch_dir / "empty.gz";

  {
    GzipOstream gz_out(gz_path);
  }  // Close() in destructor; deflate(Z_FINISH) emits header + trailer with no payload.

  // A gzip stream with zero payload still has the 10-byte header + 8-byte
  // trailer + a tiny empty deflate block; just verify we can decode it back
  // to an empty payload.
  auto const decompressed = DecompressGzip(ReadAllBytes(gz_path));
  CHECK(decompressed.empty());

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("GzipOstream: concatenates multiple Write() calls correctly",
          "[lancet][base][GzipOstream]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_gzip_ostream_multi_write");
  auto const gz_path = scratch_dir / "multi.gz";

  auto const part_one = std::string("first chunk; ");
  auto const part_two = std::string("second chunk; ");
  auto const part_three = std::string("third chunk.");

  {
    GzipOstream gz_out(gz_path);
    gz_out.Write(part_one);
    gz_out.Write(part_two);
    gz_out.Write(part_three);
  }
  auto const decompressed = DecompressGzip(ReadAllBytes(gz_path));
  // std::string ctor takes char const*; decompressed buffer is u8 (zlib byte type).
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* const decompressed_chars = reinterpret_cast<char const*>(decompressed.data());
  auto const decompressed_str = std::string{decompressed_chars, decompressed.size()};
  auto expected_payload = part_one;
  expected_payload.append(part_two);
  expected_payload.append(part_three);
  CHECK(decompressed_str == expected_payload);

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("GzipOstream: round-trips a large repetitive payload (compresses well)",
          "[lancet][base][GzipOstream]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_gzip_ostream_large");
  auto const gz_path = scratch_dir / "large.gz";

  // 1 MiB of a repeating pattern compresses heavily, so the gzip output
  // exercises the inner deflate→drain loop multiple times.
  static constexpr usize ONE_MEBIBYTE = usize{1024} * 1024;
  auto repeating_payload = std::string(ONE_MEBIBYTE, 'A');
  repeating_payload.append(ONE_MEBIBYTE, 'C');
  repeating_payload.append(ONE_MEBIBYTE, 'G');
  CHECK(RoundTripGzip(gz_path, repeating_payload) == repeating_payload);

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("GzipOstream::Close is idempotent", "[lancet][base][GzipOstream]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_gzip_ostream_close_idempotent");
  auto const gz_path = scratch_dir / "idempotent.gz";

  GzipOstream gz_out(gz_path);
  gz_out.Write("payload-bytes");
  gz_out.Close();
  // Second Close() must be a no-op (no throw, no double-trailer, no
  // double-close).
  CHECK_NOTHROW(gz_out.Close());

  auto const decompressed = DecompressGzip(ReadAllBytes(gz_path));
  // std::string ctor takes char const*; decompressed buffer is u8 (zlib byte type).
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* const decompressed_chars = reinterpret_cast<char const*>(decompressed.data());
  CHECK(std::string{decompressed_chars, decompressed.size()} == "payload-bytes");

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("GzipOstream::Write after Close throws", "[lancet][base][GzipOstream]") {
  auto const scratch_dir = MakeFreshScratchDir("lancet_gzip_ostream_write_after_close");
  auto const gz_path = scratch_dir / "post_close.gz";

  GzipOstream gz_out(gz_path);
  gz_out.Close();
  CHECK_THROWS_AS(gz_out.Write("late-bytes"), std::runtime_error);

  std::filesystem::remove_all(scratch_dir);
}

TEST_CASE("GzipOstream constructor throws when output path can't be opened",
          "[lancet][base][GzipOstream]") {
  // A path inside a non-existent directory makes ofstream::open fail.
  auto const bogus_path = std::filesystem::temp_directory_path() /
                          "lancet_gzip_ostream_no_such_dir" /
                          "subdir" /
                          "out.gz";
  std::filesystem::remove_all(bogus_path.parent_path().parent_path());

  CHECK_THROWS_AS(GzipOstream(bogus_path), std::runtime_error);
}

TEST_CASE("GzipOstream sink ctor round-trips a payload through an in-memory stringstream",
          "[lancet][base][GzipOstream]") {
  // The sink ctor lets tests inspect the gzip wire format without touching
  // the filesystem. Caller owns the stringstream — Close() flushes the gzip
  // trailer through but does NOT close the sink itself, so we can read back
  // the bytes after the GzipOstream goes out of scope.
  auto const payload = std::string("sink-ctor round-trip — no temp file needed.\n");

  std::stringstream sink(std::ios::in | std::ios::out | std::ios::binary);
  {
    GzipOstream gz_out(sink);
    gz_out.Write(payload);
  }  // ~GzipOstream flushes via mSink->flush() in this branch.

  // Pull the gzip bytes out of the sink and decompress for comparison.
  auto const gz_str = sink.str();
  std::vector<u8> gz_bytes(gz_str.size(), u8{0});
  // string::data() returns char const*; our buffer is u8 (zlib's `Bytef`) —
  // layout-compatible bytes.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* const src_bytes = reinterpret_cast<u8 const*>(gz_str.data());
  for (usize idx = 0; idx < gz_str.size(); ++idx) gz_bytes[idx] = src_bytes[idx];

  auto const decompressed = lancet::tests::DecompressGzip(gz_bytes);
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  auto const* const decompressed_chars = reinterpret_cast<char const*>(decompressed.data());
  CHECK(std::string(decompressed_chars, decompressed.size()) == payload);
}

TEST_CASE("GzipOstream sink ctor leaves the caller's sink open after Close",
          "[lancet][base][GzipOstream]") {
  // Close() in sink-ctor mode must NOT close the underlying sink — the
  // caller owns its lifetime. This guards against the sink being prematurely
  // marked closed (e.g., a regression that reused the path-ctor close path).
  std::stringstream sink(std::ios::in | std::ios::out | std::ios::binary);
  {
    GzipOstream gz_out(sink);
    gz_out.Write("payload-bytes");
    gz_out.Close();
  }
  // good() iff no fail/bad/eof bits; the sink should still accept further
  // writes if the caller wanted to append.
  CHECK(sink.good());
  CHECK_FALSE(sink.str().empty());
}
