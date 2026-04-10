// ============================================================================
// AnalyzeLongdustScores — Single-Threaded Longdust Score Analysis Tool
//
// Reads the scored BED output from ScoreRegionLongdust (potentially 500M+ rows)
// and generates two analysis files:
//   1. longdust_analysis_by_category.txt    - Median/mean/max scores per category
//   2. longdust_analysis_by_subcategory.txt - Median scores per subcategory
//
// Architecture: single-pass streaming with O(1) per-line accumulation.
//   • Reads gzipped input line-by-line via gzgets
//   • Accumulates scores into flat_hash_map accumulators (O(1) insert)
//   • Converts to btree_map only at report time for sorted output
//   • Computes exact median via nth_element at report time
//
// Input format (9-column BED from ScoreRegionLongdust):
//   #chrom  start  end  source  name  region  scale  region_length  score
//
// Usage: AnalyzeLongdustScores <scored.bed.gz> <output_dir>
// ============================================================================

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include <cmath>

extern "C" {
#include "zlib.h"
}

#include "lancet/base/timer.h"
#include "lancet/base/types.h"

#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_split.h"
#include "absl/strings/strip.h"
#include "absl/time/time.h"
#include "spdlog/fmt/bundled/format.h"

namespace {

// ── Constants ───────────────────────────────────────────────────────────────

/// Column indices in the 9-column BED format.
constexpr usize COL_SOURCE = 3;
constexpr usize COL_NAME = 4;
constexpr usize COL_SCALE = 6;
constexpr usize COL_SCORE = 8;
constexpr usize MIN_COLUMNS = 9;

/// Expected row count for pre-allocation (balanced: ≤5000/subcat × ~80 subcats × 4 scales).
constexpr usize EXPECTED_TOTAL_ROWS = 2'000'000;

/// Known dataset shape for pre-allocation.
constexpr usize NUM_SCALES = 4;
constexpr usize EXPECTED_CATEGORIES = 16;
constexpr usize EXPECTED_SUBCATS_PER = 64;

// ── Data Structures ─────────────────────────────────────────────────────────

/// Accumulated scores for one (group, scale) combination.
/// Collects all scores for exact median via nth_element at report time.
struct ScoreAccumulator {
  std::vector<double> mScores;

  explicit ScoreAccumulator(usize reserve_hint = 0) { mScores.reserve(reserve_hint); }

  void Add(double score) { mScores.push_back(score); }

  [[nodiscard]] auto Count() const -> usize { return mScores.size(); }

  [[nodiscard]] auto Median() -> double {
    if (mScores.empty()) return 0.0;
    auto const mid = mScores.size() / 2;
    std::nth_element(mScores.begin(), mScores.begin() + static_cast<i64>(mid), mScores.end());
    return mScores[mid];
  }

  [[nodiscard]] auto Mean() const -> double {
    if (mScores.empty()) return 0.0;

    double const sum = std::accumulate(mScores.cbegin(), mScores.cend(), 0.0);
    return sum / static_cast<double>(mScores.size());
  }

  [[nodiscard]] auto Max() -> double {
    if (mScores.empty()) return 0.0;
    return *std::max_element(mScores.begin(), mScores.end());
  }

  [[nodiscard]] auto CountAbove(double threshold) const -> usize {
    return static_cast<usize>(std::count_if(mScores.begin(), mScores.end(),
                                            [threshold](double scr) { return scr > threshold; }));
  }
};

/// Score data for one group (category or subcategory) across all scales.
struct GroupScores {
  absl::flat_hash_map<i64, ScoreAccumulator> mByScale;

  GroupScores() { mByScale.reserve(NUM_SCALES); }
};

/// All accumulated data for the analysis.
struct Accumulators {
  absl::flat_hash_map<std::string, GroupScores> mCategory;
  absl::flat_hash_map<std::string, absl::flat_hash_map<std::string, GroupScores>> mSubcategory;
  absl::flat_hash_map<i64, ScoreAccumulator> mGlobal;
  usize mRows = 0;

  void Reserve() {
    mCategory.reserve(EXPECTED_CATEGORIES);
    mSubcategory.reserve(EXPECTED_CATEGORIES);
    mGlobal.reserve(NUM_SCALES);
  }
};

// ── Report Writers ──────────────────────────────────────────────────────────

auto FormatScaleHeader(std::vector<i64> const& scales) -> std::string {
  std::string header;
  for (auto const scale : scales) {
    header += fmt::format("{:>9}", fmt::format("{}bp", scale * 2));
  }
  return header;
}

void WriteCategoryAnalysis(std::filesystem::path const& path,
                           absl::btree_map<std::string, GroupScores>& category_data,
                           absl::flat_hash_map<i64, ScoreAccumulator>& global_data,
                           std::vector<i64> const& scales, usize total_rows) {
  std::ofstream out(path);
  auto const scale_hdr = FormatScaleHeader(scales);

  out << fmt::format("Longdust Score Analysis by Category ({} total rows)\n", total_rows);
  out << std::string(120, '=') << "\n\n";

  // Category median table
  out << fmt::format("{:<40}{:>12}  {}\n", "Category", "N", scale_hdr);
  out << std::string(120, '-') << "\n";

  for (auto& [category, group] : category_data) {
    auto const count = group.mByScale.empty() ? 0UL : group.mByScale.begin()->second.Count();
    out << fmt::format("{:<40}{:>12}  ", category, count);
    for (auto const scale : scales) {
      out << fmt::format("{:>9.4f}", group.mByScale[scale].Median());
    }
    out << "\n";
  }

  // Global distribution
  out << "\n\nScore Distribution by Scale (across all regions)\n";
  out << std::string(120, '=') << "\n";
  for (auto const scale : scales) {
    auto& accum = global_data[scale];
    auto const count = accum.Count();
    if (count == 0) continue;

    auto const nonzero = accum.CountAbove(0.0);
    auto const gt01 = accum.CountAbove(0.1);
    auto const gt06 = accum.CountAbove(0.6);
    auto const gt10 = accum.CountAbove(1.0);
    auto const pct_nz =
        count > 0 ? 100.0 * static_cast<double>(nonzero) / static_cast<double>(count) : 0.0;
    out << fmt::format("  {:>5}bp: n={:>12}  mean={:>7.4f}  median={:>7.4f}  max={:>7.4f}"
                       "  nonzero={:>12} ({:>3.0f}%)  >0.1={:>10}  >0.6={:>10}  >1.0={:>10}\n",
                       scale * 2, count, accum.Mean(), accum.Median(), accum.Max(), nonzero, pct_nz,
                       gt01, gt06, gt10);
  }

  fmt::print(stderr, "Wrote: {}\n", path.string());
}

void WriteSubcategoryAnalysis(
    std::filesystem::path const& path,
    absl::btree_map<std::string, absl::btree_map<std::string, GroupScores>>& subcat_data,
    usize total_rows) {
  // Collect all scales from the data
  absl::flat_hash_set<i64> scale_set;
  for (auto& [cat_key, subcategories] : subcat_data) {
    for (auto& [sub_key, group] : subcategories) {
      for (auto& [scale, unused_acc] : group.mByScale) {
        scale_set.insert(scale);
      }
    }
  }
  std::vector<i64> scales(scale_set.begin(), scale_set.end());
  std::ranges::sort(scales);

  std::ofstream out(path);
  auto const scale_hdr = FormatScaleHeader(scales);

  out << fmt::format("Longdust Score Analysis by Subcategory ({} total rows)\n", total_rows);
  out << std::string(140, '=') << "\n\n";

  for (auto& [category, subcategories] : subcat_data) {
    out << fmt::format("\n── {} ({} subcategories) ──\n", category, subcategories.size());
    out << fmt::format("  {:<55}{:>10}  {}\n", "Subcategory", "N", scale_hdr);
    out << "  " << std::string(136, '-') << "\n";

    // Sort by median at largest scale (descending)
    auto const largest_scale = scales.back();
    std::vector<std::pair<double, std::string>> sorted;
    for (auto& [subcat, group] : subcategories) {
      sorted.emplace_back(-group.mByScale[largest_scale].Median(), subcat);
    }
    std::ranges::sort(sorted);

    for (auto& [neg_median, subcat] : sorted) {
      auto& group = subcategories[subcat];
      auto const count = group.mByScale.empty() ? 0UL : group.mByScale.begin()->second.Count();
      auto const display = subcat.substr(0, 55);

      out << fmt::format("  {:<55}{:>10}  ", display, count);
      for (auto const scale : scales) {
        out << fmt::format("{:>9.4f}", group.mByScale[scale].Median());
      }
      out << "\n";
    }
  }

  fmt::print(stderr, "Wrote: {}\n", path.string());
}

}  // namespace

// ── Main ────────────────────────────────────────────────────────────────────

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
auto main(int argc, char** argv) -> int {
  if (argc < 3) {
    fmt::print(stderr, "Usage: AnalyzeLongdustScores <scored.bed.gz> <output_dir>\n\n"
                       "Reads a scored BED from ScoreRegionLongdust and generates:\n"
                       "  <output_dir>/longdust_analysis_by_category.txt\n"
                       "  <output_dir>/longdust_analysis_by_subcategory.txt\n");
    return 1;
  }

  std::filesystem::path const input_path = argv[1];
  std::filesystem::path const output_dir = argv[2];
  std::filesystem::create_directories(output_dir);

  fmt::print(stderr,
             "=== AnalyzeLongdustScores ===\n"
             "Input:  {}\nOutput: {}\n\n",
             input_path.string(), output_dir.string());

  // ── Step 1: Open input ──────────────────────────────────────────────────
  gzFile gzfp = gzopen(input_path.c_str(), "r");
  if (!gzfp) {
    fmt::print(stderr, "Error: cannot open {}\n", input_path.string());
    return 1;
  }

  // ── Step 2: Stream and accumulate ───────────────────────────────────────
  Timer timer;
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
  char buf[65'536];
  Accumulators accum;
  accum.Reserve();

  fmt::print(stderr, "Streaming input...\n");
  while (gzgets(gzfp, buf, sizeof(buf)) != nullptr) {
    std::string line(buf);
    absl::StripTrailingAsciiWhitespace(&line);
    if (line.empty() || line[0] == '#' || line.starts_with("chrom")) continue;

    std::vector<std::string_view> const fields = absl::StrSplit(line, '\t');
    if (fields.size() < MIN_COLUMNS) continue;

    auto const source = std::string(fields[COL_SOURCE]);
    auto const fname = std::string(fields[COL_NAME]);
    auto const scale = std::stol(std::string(fields[COL_SCALE]));
    auto const score = std::stod(std::string(fields[COL_SCORE]));

    accum.mCategory[source].mByScale[scale].Add(score);
    accum.mSubcategory[source][fname].mByScale[scale].Add(score);
    accum.mGlobal[scale].Add(score);
    accum.mRows++;

    if (accum.mRows % 10'000'000 == 0) {
      auto const elapsed = absl::FormatDuration(absl::Trunc(timer.Runtime(), absl::Seconds(1)));
      fmt::print(stderr, "  {}M rows processed | {}\n", accum.mRows / 1'000'000, elapsed);
    }
  }
  gzclose(gzfp);

  auto const stream_elapsed = absl::FormatDuration(absl::Trunc(timer.Runtime(), absl::Seconds(1)));
  fmt::print(stderr, "Streamed {} rows ({})\n\n", accum.mRows, stream_elapsed);

  // ── Step 3: Determine scales and convert to btree_map for sorted output ─
  std::vector<i64> scales;
  for (auto const& [scale, unused_acc] : accum.mGlobal) {
    scales.push_back(scale);
  }
  std::ranges::sort(scales);

  absl::btree_map<std::string, GroupScores> sorted_categories(accum.mCategory.begin(),
                                                              accum.mCategory.end());

  absl::btree_map<std::string, absl::btree_map<std::string, GroupScores>> sorted_subcategories;
  for (auto& [cat_key, subcats] : accum.mSubcategory) {
    sorted_subcategories[cat_key] =
        absl::btree_map<std::string, GroupScores>(subcats.begin(), subcats.end());
  }

  // ── Step 4: Write analysis reports ──────────────────────────────────────
  fmt::print(stderr, "Writing analysis reports...\n");
  WriteCategoryAnalysis(output_dir / "longdust_analysis_by_category.txt", sorted_categories,
                        accum.mGlobal, scales, accum.mRows);

  WriteSubcategoryAnalysis(output_dir / "longdust_analysis_by_subcategory.txt",
                           sorted_subcategories, accum.mRows);

  fmt::print(stderr, "\n=== Done: analyzed {} rows ===\n", accum.mRows);
  return 0;
}
