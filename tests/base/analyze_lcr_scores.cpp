// ============================================================================
// AnalyzeLCRScores — Single-Threaded LCR Score Analysis Tool
//
// Reads the scored BED output from ScoreRegionLCR (potentially 500M+ rows)
// and generates two analysis files:
//   1. lcr_analysis_by_category.txt    - Median/mean/max scores per category
//   2. lcr_analysis_by_subcategory.txt - Median scores per subcategory
//
// Architecture: single-pass streaming with O(1) per-line accumulation.
//   • Reads gzipped input line-by-line via gzgets
//   • Accumulates scores into flat_hash_map accumulators (O(1) insert)
//   • Converts to btree_map only at report time for sorted output
//   • Computes exact median via nth_element at report time
//
// Input format (9-column BED from ScoreRegionLCR):
//   #chrom  start  end  source  name  region  scale  region_length  score
//
// Usage: AnalyzeLCRScores <scored.bed.gz> <output_dir>
// ============================================================================

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

extern "C" {
#include "zlib.h"
}

#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_split.h"
#include "absl/strings/strip.h"
#include "absl/time/time.h"
#include "spdlog/fmt/bundled/format.h"

#include "lancet/base/timer.h"
#include "lancet/base/types.h"

namespace {

// ── Constants ───────────────────────────────────────────────────────────────

/// Column indices in the 9-column BED format.
static constexpr usize COL_SOURCE = 3;
static constexpr usize COL_NAME = 4;
static constexpr usize COL_SCALE = 6;
static constexpr usize COL_SCORE = 8;
static constexpr usize MIN_COLUMNS = 9;

/// Expected row count for pre-allocation (balanced: ≤5000/subcat × ~80 subcats × 4 scales).
static constexpr usize EXPECTED_TOTAL_ROWS = 2'000'000;

/// Known dataset shape for pre-allocation.
static constexpr usize NUM_SCALES = 4;
static constexpr usize EXPECTED_CATEGORIES = 16;
static constexpr usize EXPECTED_SUBCATS_PER = 64;

// ── Data Structures ─────────────────────────────────────────────────────────

/// Accumulated scores for one (group, scale) combination.
/// Collects all scores for exact median via nth_element at report time.
struct ScoreAccumulator {
  std::vector<double> scores;

  explicit ScoreAccumulator(usize reserve_hint = 0) { scores.reserve(reserve_hint); }

  void Add(double score) { scores.push_back(score); }

  [[nodiscard]] auto Count() const -> usize { return scores.size(); }

  [[nodiscard]] auto Median() -> double {
    if (scores.empty()) return 0.0;
    const auto mid = scores.size() / 2;
    std::nth_element(scores.begin(), scores.begin() + static_cast<i64>(mid), scores.end());
    return scores[mid];
  }

  [[nodiscard]] auto Mean() const -> double {
    if (scores.empty()) return 0.0;
    double sum = 0;
    for (const auto val : scores) sum += val;
    return sum / static_cast<double>(scores.size());
  }

  [[nodiscard]] auto Max() -> double {
    if (scores.empty()) return 0.0;
    return *std::max_element(scores.begin(), scores.end());
  }

  [[nodiscard]] auto CountAbove(double threshold) const -> usize {
    return static_cast<usize>(
        std::count_if(scores.begin(), scores.end(), [threshold](double s) { return s > threshold; }));
  }
};

/// Score data for one group (category or subcategory) across all scales.
struct GroupScores {
  absl::flat_hash_map<i64, ScoreAccumulator> by_scale;

  GroupScores() { by_scale.reserve(NUM_SCALES); }
};

/// All accumulated data for the analysis.
struct Accumulators {
  absl::flat_hash_map<std::string, GroupScores> category;
  absl::flat_hash_map<std::string, absl::flat_hash_map<std::string, GroupScores>> subcategory;
  absl::flat_hash_map<i64, ScoreAccumulator> global;
  usize rows = 0;

  void Reserve() {
    category.reserve(EXPECTED_CATEGORIES);
    subcategory.reserve(EXPECTED_CATEGORIES);
    global.reserve(NUM_SCALES);
  }
};

// ── Report Writers ──────────────────────────────────────────────────────────

auto FormatScaleHeader(const std::vector<i64>& scales) -> std::string {
  std::string hdr;
  for (const auto scale : scales) {
    hdr += fmt::format("{:>9}", fmt::format("{}bp", scale * 2));
  }
  return hdr;
}

void WriteCategoryAnalysis(const std::filesystem::path& path,
                           absl::btree_map<std::string, GroupScores>& category_data,
                           absl::flat_hash_map<i64, ScoreAccumulator>& global_data,
                           const std::vector<i64>& scales,
                           usize total_rows) {
  std::ofstream out(path);
  const auto scale_hdr = FormatScaleHeader(scales);

  out << fmt::format("LCR Score Analysis by Category ({} total rows)\n", total_rows);
  out << std::string(120, '=') << "\n\n";

  // Category median table
  out << fmt::format("{:<40}{:>12}  {}\n", "Category", "N", scale_hdr);
  out << std::string(120, '-') << "\n";

  for (auto& [category, group] : category_data) {
    const auto n = group.by_scale.empty() ? 0UL : group.by_scale.begin()->second.Count();
    out << fmt::format("{:<40}{:>12}  ", category, n);
    for (const auto scale : scales) {
      out << fmt::format("{:>9.4f}", group.by_scale[scale].Median());
    }
    out << "\n";
  }

  // Global distribution
  out << "\n\nScore Distribution by Scale (across all regions)\n";
  out << std::string(120, '=') << "\n";
  for (const auto scale : scales) {
    auto& acc = global_data[scale];
    const auto n = acc.Count();
    if (n == 0) continue;
    const auto nonzero = acc.CountAbove(0.0);
    const auto gt01 = acc.CountAbove(0.1);
    const auto gt06 = acc.CountAbove(0.6);
    const auto gt10 = acc.CountAbove(1.0);
    const auto pct_nz = n > 0 ? 100.0 * static_cast<double>(nonzero) / static_cast<double>(n) : 0.0;
    out << fmt::format("  {:>5}bp: n={:>12}  mean={:>7.4f}  median={:>7.4f}  max={:>7.4f}"
                       "  nonzero={:>12} ({:>3.0f}%)  >0.1={:>10}  >0.6={:>10}  >1.0={:>10}\n",
                       scale * 2, n, acc.Mean(), acc.Median(), acc.Max(),
                       nonzero, pct_nz, gt01, gt06, gt10);
  }

  fmt::print(stderr, "Wrote: {}\n", path.string());
}

void WriteSubcategoryAnalysis(
    const std::filesystem::path& path,
    absl::btree_map<std::string, absl::btree_map<std::string, GroupScores>>& subcat_data,
    usize total_rows) {
  // Collect all scales from the data
  absl::flat_hash_set<i64> scale_set;
  for (auto& [_, subcategories] : subcat_data) {
    for (auto& [__, group] : subcategories) {
      for (auto& [scale, ___] : group.by_scale) {
        scale_set.insert(scale);
      }
    }
  }
  std::vector<i64> scales(scale_set.begin(), scale_set.end());
  std::ranges::sort(scales);

  std::ofstream out(path);
  const auto scale_hdr = FormatScaleHeader(scales);

  out << fmt::format("LCR Score Analysis by Subcategory ({} total rows)\n", total_rows);
  out << std::string(140, '=') << "\n\n";

  for (auto& [category, subcategories] : subcat_data) {
    out << fmt::format("\n── {} ({} subcategories) ──\n", category, subcategories.size());
    out << fmt::format("  {:<55}{:>10}  {}\n", "Subcategory", "N", scale_hdr);
    out << "  " << std::string(136, '-') << "\n";

    // Sort by median at largest scale (descending)
    const auto largest_scale = scales.back();
    std::vector<std::pair<double, std::string>> sorted;
    for (auto& [subcat, group] : subcategories) {
      sorted.emplace_back(-group.by_scale[largest_scale].Median(), subcat);
    }
    std::ranges::sort(sorted);

    for (auto& [_, subcat] : sorted) {
      auto& group = subcategories[subcat];
      const auto n = group.by_scale.empty() ? 0UL : group.by_scale.begin()->second.Count();
      const auto display = subcat.substr(0, 55);

      out << fmt::format("  {:<55}{:>10}  ", display, n);
      for (const auto scale : scales) {
        out << fmt::format("{:>9.4f}", group.by_scale[scale].Median());
      }
      out << "\n";
    }
  }

  fmt::print(stderr, "Wrote: {}\n", path.string());
}

}  // namespace

// ── Main ────────────────────────────────────────────────────────────────────

auto main(int argc, char** argv) -> int {
  if (argc < 3) {
    fmt::print(stderr,
               "Usage: AnalyzeLCRScores <scored.bed.gz> <output_dir>\n\n"
               "Reads a scored BED from ScoreRegionLCR and generates:\n"
               "  <output_dir>/lcr_analysis_by_category.txt\n"
               "  <output_dir>/lcr_analysis_by_subcategory.txt\n");
    return 1;
  }

  const std::filesystem::path input_path = argv[1];
  const std::filesystem::path output_dir = argv[2];
  std::filesystem::create_directories(output_dir);

  fmt::print(stderr,
             "=== AnalyzeLCRScores ===\n"
             "Input:  {}\nOutput: {}\n\n",
             input_path.string(), output_dir.string());

  // ── Step 1: Open input ──────────────────────────────────────────────────
  gzFile gz = gzopen(input_path.c_str(), "r");
  if (!gz) {
    fmt::print(stderr, "Error: cannot open {}\n", input_path.string());
    return 1;
  }

  // ── Step 2: Stream and accumulate ───────────────────────────────────────
  Timer timer;
  char buf[65536];
  Accumulators acc;
  acc.Reserve();

  fmt::print(stderr, "Streaming input...\n");
  while (gzgets(gz, buf, sizeof(buf)) != nullptr) {
    std::string line(buf);
    absl::StripTrailingAsciiWhitespace(&line);
    if (line.empty() || line[0] == '#' || line.starts_with("chrom")) continue;

    const std::vector<std::string_view> fields = absl::StrSplit(line, '\t');
    if (fields.size() < MIN_COLUMNS) continue;

    const auto source = std::string(fields[COL_SOURCE]);
    const auto name = std::string(fields[COL_NAME]);
    const auto scale = std::stol(std::string(fields[COL_SCALE]));
    const auto score = std::stod(std::string(fields[COL_SCORE]));

    acc.category[source].by_scale[scale].Add(score);
    acc.subcategory[source][name].by_scale[scale].Add(score);
    acc.global[scale].Add(score);
    acc.rows++;

    if (acc.rows % 10'000'000 == 0) {
      const auto elapsed = absl::FormatDuration(absl::Trunc(timer.Runtime(), absl::Seconds(1)));
      fmt::print(stderr, "  {}M rows processed | {}\n", acc.rows / 1'000'000, elapsed);
    }
  }
  gzclose(gz);

  const auto stream_elapsed = absl::FormatDuration(absl::Trunc(timer.Runtime(), absl::Seconds(1)));
  fmt::print(stderr, "Streamed {} rows ({})\n\n", acc.rows, stream_elapsed);

  // ── Step 3: Determine scales and convert to btree_map for sorted output ─
  std::vector<i64> scales;
  for (const auto& [scale, _] : acc.global) {
    scales.push_back(scale);
  }
  std::ranges::sort(scales);

  absl::btree_map<std::string, GroupScores> sorted_categories(
      acc.category.begin(), acc.category.end());

  absl::btree_map<std::string, absl::btree_map<std::string, GroupScores>> sorted_subcategories;
  for (auto& [cat, subcats] : acc.subcategory) {
    sorted_subcategories[cat] = absl::btree_map<std::string, GroupScores>(
        subcats.begin(), subcats.end());
  }

  // ── Step 4: Write analysis reports ──────────────────────────────────────
  fmt::print(stderr, "Writing analysis reports...\n");
  WriteCategoryAnalysis(output_dir / "lcr_analysis_by_category.txt",
                        sorted_categories, acc.global, scales, acc.rows);

  WriteSubcategoryAnalysis(output_dir / "lcr_analysis_by_subcategory.txt",
                           sorted_subcategories, acc.rows);

  fmt::print(stderr, "\n=== Done: analyzed {} rows ===\n", acc.rows);
  return 0;
}
