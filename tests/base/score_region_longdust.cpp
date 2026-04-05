// ============================================================================
// ScoreRegionLongdust — Comprehensive Longdust Score Analysis & Test Data Generator
//
// Reads annotation BED files in parallel, expands each to 4 scoring windows
// (one per scale), sorts all windows globally by genome position, scores them
// in parallel, and writes two output files:
//   1. BGZF-compressed scored BED (balanced: ≤5000 regions per subcategory)
//   2. Test calibration TSV (deterministically sampled ~8 per subcategory)
//
// Architecture (modeled after pipeline_runner.cpp):
//   • All chromosome sequences cached in memory for O(1) random access
//   • BED annotation files loaded in parallel (especially GIAB gzipped files)
//   • Each annotation expanded to 8 ScoringWindows, then ALL windows sorted
//     globally by (chrom, window_start, window_end, source, name, scale)
//   • moodycamel::ConcurrentQueue for work distribution and result collection
//   • Ordered output flushing (process in any order, write in genome order)
//   • Deterministic test sampling pre-computed per subcategory
//
// Output BED columns (0-based half-open):
//   chrom  start  end  source  name  region  scale  region_length  score
//
// Usage:
//   ScoreRegionLongdust <ref.fa.gz> <data_dir> <output.bed.gz> <test.tsv> <threads>
// ============================================================================

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <stop_token>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

extern "C" {
#include "zlib.h"
}

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/fixed_array.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_split.h"
#include "absl/strings/strip.h"
#include "absl/synchronization/mutex.h"
#include "absl/time/time.h"
#include "concurrentqueue.h"
#include "spdlog/fmt/bundled/format.h"

#include "lancet/base/longdust_scorer.h"
#include "lancet/base/timer.h"
#include "lancet/base/types.h"
#include "lancet/hts/bgzf_ostream.h"
#include "lancet/hts/reference.h"

namespace {

// ── Constants ───────────────────────────────────────────────────────────────

static constexpr std::array<i64, 4> SCALES = {50, 100, 500, 1000};
static constexpr usize NUM_SCALES = SCALES.size();
static constexpr usize BALANCED_SAMPLES_PER_SUBCAT = 5000;
static constexpr usize TEST_SAMPLES_PER_SUBCAT = 8;

static constexpr std::array<std::string_view, 12> GIAB_SKIP_PREFIXES = {
    "AllHomopolymers",      "AllTandemRepeats",        "allTandemRepeats",
    "notinAllHomopolymers", "notinAllTandemRepeats",   "notinallTandemRepeats",
    "notinsatellites",      "alldifficultregions",     "notinSegDups",
    "AllAutosomes",         "chrX_",                   "chrY_",
};

// ── Data Structures ─────────────────────────────────────────────────────────

/// A loaded BED annotation with source-specific category/subcategory.
struct BedAnnotation {
  std::string chrom;
  i64 start;           // 0-based inclusive
  i64 end;             // 0-based exclusive
  std::string source;  // e.g. "RepeatMasker", "GIAB/LowComplexity"
  std::string name;    // subcategory (source-specific column)
  u8 chrom_idx = 0;    // FASTA chromosome order index
};

/// A single scoring window derived from one BedAnnotation at one scale.
/// This is the atomic unit of work, sorting, and output.
struct ScoringWindow {
  usize parent_idx;    // index into annotations[] (for source/name/chrom lookup)
  u8 chrom_idx;        // FASTA order, duplicated for fast sorting
  i64 window_start;    // 0-based inclusive
  i64 window_end;      // 0-based exclusive
  i64 scale;
  i64 region_length;   // original annotation span (end - start)
};

/// Result from scoring a single window.
struct ScoredWindow {
  usize window_idx;    // index into windows[]
  double score;
};

// ── Utility ─────────────────────────────────────────────────────────────────

auto IsValidChrom(std::string_view chrom) -> bool {
  return chrom.find("Un") == std::string_view::npos && chrom != "chrM";
}

auto ShouldSkipGiab(std::string_view name) -> bool {
  return std::ranges::any_of(GIAB_SKIP_PREFIXES,
                             [&](auto pfx) { return name.starts_with(pfx); });
}

// ── Reference Caching ───────────────────────────────────────────────────────

struct CachedReference {
  absl::flat_hash_map<std::string, std::string> sequences;
  absl::flat_hash_map<std::string, i64> lengths;
  absl::flat_hash_map<std::string, u8> chrom_order;

  static auto Build(const lancet::hts::Reference& ref) -> CachedReference {
    CachedReference cache;
    const auto chroms = ref.ListChroms();
    fmt::print(stderr, "Caching {} chromosome sequences...\n", chroms.size());
    for (const auto& chrom : chroms) {
      const auto nm = chrom.Name();
      const auto len = static_cast<i64>(chrom.Length());
      fmt::print(stderr, "  {} ({} bp)\n", nm, len);
      cache.sequences.emplace(nm, std::string(ref.MakeRegion(nm.c_str()).SeqView()));
      cache.lengths.emplace(nm, len);
      cache.chrom_order.emplace(nm, static_cast<u8>(chrom.Index()));
    }
    fmt::print(stderr, "Done caching reference.\n\n");
    return cache;
  }
};

// ── BED Parsing (source-specific) ───────────────────────────────────────────

/// Parse a plain-text BED, assigning source/name per the caller's spec.
/// subcat_column: 0-indexed column index for per-row subcategory (-1 = use fixed_name).
auto ParsePlainBed(const std::filesystem::path& filepath,
                   const CachedReference& cache,
                   const std::string& source,
                   const std::string& fixed_name,
                   int subcat_column = -1) -> std::vector<BedAnnotation> {
  std::vector<BedAnnotation> result;
  std::ifstream file(filepath);
  if (!file.is_open()) {
    fmt::print(stderr, "  Warning: cannot open {}\n", filepath.string());
    return result;
  }

  std::string line;
  while (std::getline(file, line)) {
    absl::StripTrailingAsciiWhitespace(&line);
    if (line.empty() || line[0] == '#' || line[0] == 't') continue;

    const std::vector<std::string_view> fields = absl::StrSplit(line, '\t');
    if (fields.size() < 3) continue;
    const auto chrom = std::string(fields[0]);
    if (!IsValidChrom(chrom)) continue;

    const auto it = cache.chrom_order.find(chrom);
    if (it == cache.chrom_order.end()) continue;

    std::string name = fixed_name;
    if (subcat_column >= 0 && static_cast<usize>(subcat_column) < fields.size()) {
      name = std::string(fields[subcat_column]);
    }

    result.push_back({chrom, std::stol(std::string(fields[1])),
                      std::stol(std::string(fields[2])),
                      source, std::move(name), it->second});
  }
  return result;
}

/// Parse a gzipped BED (for GIAB .bed.gz files).
auto ParseGzippedBed(const std::filesystem::path& filepath,
                     const CachedReference& cache,
                     const std::string& source,
                     const std::string& fixed_name) -> std::vector<BedAnnotation> {
  std::vector<BedAnnotation> result;
  gzFile gz = gzopen(filepath.c_str(), "r");
  if (!gz) {
    fmt::print(stderr, "  Warning: cannot open {}\n", filepath.string());
    return result;
  }

  char buf[4096];
  while (gzgets(gz, buf, sizeof(buf)) != nullptr) {
    std::string line(buf);
    absl::StripTrailingAsciiWhitespace(&line);
    if (line.empty() || line[0] == '#' || line[0] == 't') continue;

    const std::vector<std::string_view> fields = absl::StrSplit(line, '\t');
    if (fields.size() < 3) continue;
    const auto chrom = std::string(fields[0]);
    if (!IsValidChrom(chrom)) continue;

    const auto it = cache.chrom_order.find(chrom);
    if (it == cache.chrom_order.end()) continue;

    result.push_back({chrom, std::stol(std::string(fields[1])),
                      std::stol(std::string(fields[2])),
                      source, fixed_name, it->second});
  }
  gzclose(gz);
  return result;
}

// ── Source-Specific Loaders ─────────────────────────────────────────────────

auto LoadGiab(const std::filesystem::path& data_dir,
              const CachedReference& cache) -> std::vector<BedAnnotation> {
  std::vector<BedAnnotation> all;
  const auto giab_dir = data_dir / "chm13_giab_genome_stratifications";
  if (!std::filesystem::exists(giab_dir)) return all;

  struct GiabTarget {
    std::filesystem::path path;
    std::string source;
    std::string strat_name;
  };

  std::vector<GiabTarget> targets;
  for (const auto& entry : std::filesystem::recursive_directory_iterator(giab_dir)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".gz" ||
        entry.path().stem().extension() != ".bed") {
      continue;
    }
    auto strat = entry.path().stem().stem().string();
    const auto pfx = strat.find("CHM13v2.0_");
    if (pfx != std::string::npos) strat = strat.substr(pfx + 10);
    if (ShouldSkipGiab(strat)) continue;

    targets.push_back({entry.path(),
                       absl::StrCat("GIAB/", entry.path().parent_path().filename().string()),
                       std::move(strat)});
  }
  std::ranges::sort(targets, [](const auto& a, const auto& b) { return a.path < b.path; });
  fmt::print(stderr, "  GIAB: {} primary BED files to load in parallel\n", targets.size());

  absl::Mutex mtx;
  {
    std::vector<std::jthread> loaders;
    loaders.reserve(targets.size());
    for (const auto& tgt : targets) {
      loaders.emplace_back([&tgt, &cache, &all, &mtx] {
        auto regions = ParseGzippedBed(tgt.path, cache, tgt.source, tgt.strat_name);
        fmt::print(stderr, "    {}/{} → {} regions\n", tgt.source, tgt.strat_name, regions.size());
        absl::MutexLock lock(mtx);
        all.insert(all.end(), std::make_move_iterator(regions.begin()),
                   std::make_move_iterator(regions.end()));
      });
    }
  }

  fmt::print(stderr, "  GIAB total: {} regions\n", all.size());
  return all;
}

auto LoadRepeatMasker(const std::filesystem::path& data_dir,
                      const CachedReference& cache) -> std::vector<BedAnnotation> {
  // RepeatMasker: source="RepeatMasker", name=column 7 (0-indexed col 6, repeat class)
  return ParsePlainBed(data_dir / "chm13_repeat_annotations" /
                           "chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed",
                       cache, "RepeatMasker", "", 6);
}

auto LoadCenSat(const std::filesystem::path& data_dir,
                const CachedReference& cache) -> std::vector<BedAnnotation> {
  // CenSat: source="CenSat", name="CenSat" (fixed)
  return ParsePlainBed(data_dir / "chm13_repeat_annotations" / "chm13v2.0_censat_v2.1.bed",
                       cache, "CenSat", "CenSat");
}

auto LoadTelomere(const std::filesystem::path& data_dir,
                  const CachedReference& cache) -> std::vector<BedAnnotation> {
  // Telomere: 3-column BED, source="Telomere", name="{chrom}_{arm}_telomere"
  std::vector<BedAnnotation> result;
  const auto telo_path = data_dir / "chm13_repeat_annotations" / "chm13v2.0_telomere.bed";
  std::ifstream file(telo_path);
  std::string line;
  while (std::getline(file, line)) {
    absl::StripTrailingAsciiWhitespace(&line);
    if (line.empty()) continue;
    const std::vector<std::string_view> fields = absl::StrSplit(line, '\t');
    if (fields.size() < 3) continue;
    const auto chrom = std::string(fields[0]);
    if (!IsValidChrom(chrom)) continue;
    const auto it = cache.chrom_order.find(chrom);
    if (it == cache.chrom_order.end()) continue;

    const auto start = std::stol(std::string(fields[1]));
    const auto end = std::stol(std::string(fields[2]));
    const auto* arm = (start == 0) ? "p-arm" : "q-arm";
    result.push_back({chrom, start, end, "Telomere",
                      absl::StrCat(chrom, "_", arm, "_telomere"), it->second});
  }
  return result;
}

auto LoadNewSatellite(const std::filesystem::path& data_dir,
                      const CachedReference& cache) -> std::vector<BedAnnotation> {
  // NewSatellite: source="NewSatellite", name=column 4 (0-indexed col 3)
  return ParsePlainBed(data_dir / "chm13_repeat_annotations" /
                           "chm13v2.0_new-satellites_2022DEC.bed",
                       cache, "NewSatellite", "", 3);
}

auto LoadCompositeRepeats(const std::filesystem::path& data_dir,
                          const CachedReference& cache) -> std::vector<BedAnnotation> {
  // CompositeRepeats: source="CompositeRepeats", name=column 4 (0-indexed col 3)
  return ParsePlainBed(data_dir / "chm13_repeat_annotations" /
                           "chm13v2.0_composite-repeats_2022DEC.bed",
                       cache, "CompositeRepeats", "", 3);
}

// ── Parallel Annotation Loader ──────────────────────────────────────────────

auto LoadAllAnnotations(const std::filesystem::path& data_dir,
                        const CachedReference& cache) -> std::vector<BedAnnotation> {
  std::vector<BedAnnotation> all;
  absl::Mutex mtx;

  auto merge = [&](std::vector<BedAnnotation> batch, std::string_view label) {
    fmt::print(stderr, "  {} → {} regions\n", label, batch.size());
    absl::MutexLock lock(mtx);
    all.insert(all.end(), std::make_move_iterator(batch.begin()),
               std::make_move_iterator(batch.end()));
  };

  fmt::print(stderr, "Loading all BED annotation sources in parallel...\n");
  {
    std::vector<std::jthread> loaders;
    loaders.emplace_back([&] { merge(LoadGiab(data_dir, cache), "GIAB"); });
    loaders.emplace_back([&] { merge(LoadRepeatMasker(data_dir, cache), "RepeatMasker"); });
    loaders.emplace_back([&] { merge(LoadCenSat(data_dir, cache), "CenSat"); });
    loaders.emplace_back([&] { merge(LoadTelomere(data_dir, cache), "Telomere"); });
    loaders.emplace_back([&] { merge(LoadNewSatellite(data_dir, cache), "NewSatellite"); });
    loaders.emplace_back([&] { merge(LoadCompositeRepeats(data_dir, cache), "CompositeRepeats"); });
  }

  fmt::print(stderr, "Total annotations loaded: {}\n\n", all.size());
  return all;
}

// ── Expand Annotations to Globally Sorted Windows ───────────────────────────

/// Expand each annotation into up to NUM_SCALES windows, then sort all windows
/// globally by genome position so the output BED file is fully coordinate-sorted.
auto ExpandAndSortWindows(const std::vector<BedAnnotation>& annotations,
                          const CachedReference& cache) -> std::vector<ScoringWindow> {
  fmt::print(stderr, "Expanding {} annotations to scoring windows...\n", annotations.size());

  std::vector<ScoringWindow> windows;
  windows.reserve(annotations.size() * NUM_SCALES);

  for (usize ann_idx = 0; ann_idx < annotations.size(); ++ann_idx) {
    const auto& ann = annotations[ann_idx];
    const auto chrom_len = cache.lengths.at(ann.chrom);
    const i64 center = (ann.start + ann.end) / 2;
    const i64 region_length = ann.end - ann.start;

    for (const auto scale : SCALES) {
      const i64 ws = std::max<i64>(0, center - scale);
      const i64 we = std::min(chrom_len, center + scale);
      if (we - ws < 7) continue;

      windows.push_back({ann_idx, ann.chrom_idx, ws, we, scale, region_length});
    }
  }

  fmt::print(stderr, "Sorting {} windows by genome position...\n", windows.size());
  std::ranges::sort(windows, [](const ScoringWindow& a, const ScoringWindow& b) {
    if (a.chrom_idx != b.chrom_idx) return a.chrom_idx < b.chrom_idx;
    if (a.window_start != b.window_start) return a.window_start < b.window_start;
    if (a.window_end != b.window_end) return a.window_end < b.window_end;
    return a.scale < b.scale;
  });

  fmt::print(stderr, "Done: {} sorted scoring windows.\n\n", windows.size());
  return windows;
}

// ── Test Sample Set (deterministic, pre-computed) ───────────────────────────

/// Build a set of parent annotation indices for test TSV output.
/// Selects ~TEST_SAMPLES_PER_SUBCAT evenly-spaced annotations per (source, name).
/// All windows belonging to a sampled annotation will be included in the test file.
auto BuildTestSampleSet(const std::vector<BedAnnotation>& annotations)
    -> absl::flat_hash_set<usize> {
  // Group annotation indices by subcategory key
  absl::flat_hash_map<std::string, std::vector<usize>> groups;
  for (usize i = 0; i < annotations.size(); ++i) {
    const auto key = absl::StrCat(annotations[i].source, "\t", annotations[i].name);
    groups[key].push_back(i);
  }

  absl::flat_hash_set<usize> test_set;
  for (const auto& [key, indices] : groups) {
    const usize n = indices.size();
    const usize target = std::min(n, TEST_SAMPLES_PER_SUBCAT);
    for (usize k = 0; k < target; ++k) {
      test_set.insert(indices[k * n / target]);  // evenly spaced
    }
  }

  fmt::print(stderr, "Test sample set: {} subcategories → {} sampled annotations\n\n",
             groups.size(), test_set.size());
  return test_set;
}

// ── Worker Function ─────────────────────────────────────────────────────────

/// Input/output queue types. Work items are annotation indices.
using InputQueue = moodycamel::ConcurrentQueue<usize>;
using OutputQueue = moodycamel::ConcurrentQueue<ScoredWindow>;

/// Worker: dequeues a window index, scores the subsequence, enqueues result.
void ScoringWorker(std::stop_token stoken,
                   const moodycamel::ProducerToken& in_token,
                   std::shared_ptr<InputQueue> in_queue,
                   std::shared_ptr<OutputQueue> out_queue,
                   const std::vector<ScoringWindow>* windows,
                   const std::vector<BedAnnotation>* annotations,
                   const CachedReference* cache) {
  lancet::base::LongdustQScorer scorer(7, 10001);
  const moodycamel::ProducerToken out_token(*out_queue);
  usize win_idx = 0;

  while (true) {
    if (stoken.stop_requested()) break;
    if (!in_queue->try_dequeue_from_producer(in_token, win_idx)) continue;

    const auto& win = (*windows)[win_idx];
    const auto& ann = (*annotations)[win.parent_idx];
    const auto& chrom_seq = cache->sequences.at(ann.chrom);

    const auto len = static_cast<usize>(win.window_end - win.window_start);
    const auto subseq = std::string_view(chrom_seq).substr(static_cast<usize>(win.window_start), len);

    out_queue->enqueue(out_token, ScoredWindow{win_idx, scorer.Score(subseq)});
  }
}

// ── Output Formatting ───────────────────────────────────────────────────────

auto FormatLine(const BedAnnotation& ann, const ScoringWindow& win, double score) -> std::string {
  return fmt::format("{}\t{}\t{}\t{}\t{}\t{}:{}-{}\t{}\t{}\t{}\n",
                     ann.chrom, win.window_start, win.window_end,
                     ann.source, ann.name,
                     ann.chrom, win.window_start + 1, win.window_end,
                     win.scale, win.region_length, score);
}

// ── ETA Timer ───────────────────────────────────────────────────────────────

class EtaTimer {
 public:
  explicit EtaTimer(usize total) : total_(total) {}

  void Increment() {
    done_++;
    const double secs = absl::ToDoubleSeconds(timer_.Runtime());
    if (secs > 0) rate_ = static_cast<double>(done_) / secs;
  }

  [[nodiscard]] auto Elapsed() -> absl::Duration { return timer_.Runtime(); }

  [[nodiscard]] auto EstimatedEta() const -> absl::Duration {
    if (rate_ <= 0) return absl::InfiniteDuration();
    return absl::Seconds(static_cast<double>(total_ - done_) / rate_);
  }

  [[nodiscard]] auto RatePerSecond() const -> double { return rate_; }
  [[nodiscard]] auto PercentDone() const -> double {
    return 100.0 * static_cast<double>(done_) / static_cast<double>(total_);
  }

 private:
  usize done_ = 0;
  usize total_;
  Timer timer_;
  double rate_ = 0;
};

}  // namespace

// ── Main ────────────────────────────────────────────────────────────────────

auto main(int argc, char** argv) -> int {
  if (argc < 6) {
    fmt::print(stderr,
               "Usage: ScoreRegionLongdust <ref.fa.gz> <data_dir> <output.bed.gz> <test.tsv> <threads>\n\n"
               "Arguments:\n"
               "  ref.fa.gz      Bgzipped reference FASTA\n"
               "  data_dir       Data directory with BED annotation files\n"
               "  output.bed.gz  BGZF-compressed scored BED output (full dataset)\n"
               "  test.tsv       Test calibration TSV (sampled for regression tests)\n"
               "  threads        Number of worker threads (e.g. 96)\n");
    return 1;
  }

  const std::filesystem::path ref_path = argv[1];
  const std::filesystem::path data_dir = argv[2];
  const std::filesystem::path output_path = argv[3];
  const std::filesystem::path test_path = argv[4];
  const auto num_threads = static_cast<usize>(std::stoi(argv[5]));

  fmt::print(stderr,
             "=== ScoreRegionLongdust ===\n"
             "Reference: {}\n"
             "Data dir:  {}\n"
             "Output:    {}\n"
             "Test TSV:  {}\n"
             "Threads:   {}\n"
             "Scales:    [50, 100, 500, 1000]\n"
             "Balanced:  {} per subcategory\n\n",
             ref_path.string(), data_dir.string(), output_path.string(),
             test_path.string(), num_threads, BALANCED_SAMPLES_PER_SUBCAT);

  // ── Step 1: Cache reference genome ──────────────────────────────────────
  fmt::print(stderr, "Step 1: Loading and caching reference genome...\n");
  const lancet::hts::Reference ref(ref_path);
  auto cache = CachedReference::Build(ref);

  // ── Step 2: Load all BED annotations in parallel ────────────────────────
  fmt::print(stderr, "Step 2: Loading BED annotations...\n");
  auto all_annotations = LoadAllAnnotations(data_dir, cache);

  // ── Step 3: Balanced sampling (≤N per subcategory) ──────────────────────
  fmt::print(stderr, "Step 3: Balanced sampling (≤{} per subcategory)...\n",
             BALANCED_SAMPLES_PER_SUBCAT);
  absl::flat_hash_map<std::string, std::vector<usize>> subcat_groups;
  for (usize i = 0; i < all_annotations.size(); ++i) {
    const auto key = absl::StrCat(all_annotations[i].source, "\t", all_annotations[i].name);
    subcat_groups[key].push_back(i);
  }

  absl::flat_hash_set<usize> balanced_set;
  for (const auto& [key, indices] : subcat_groups) {
    const usize n = indices.size();
    const usize target = std::min(n, BALANCED_SAMPLES_PER_SUBCAT);
    for (usize k = 0; k < target; ++k) {
      balanced_set.insert(indices[k * n / target]);  // evenly spaced
    }
  }

  // Build the balanced annotation vector (only sampled annotations)
  std::vector<BedAnnotation> annotations;
  annotations.reserve(balanced_set.size());
  for (usize i = 0; i < all_annotations.size(); ++i) {
    if (balanced_set.contains(i)) {
      annotations.push_back(std::move(all_annotations[i]));
    }
  }
  all_annotations.clear();
  all_annotations.shrink_to_fit();

  fmt::print(stderr, "  {} subcategories → {} sampled annotations (from {} total)\n\n",
             subcat_groups.size(), annotations.size(), balanced_set.size());

  // ── Step 4: Expand to scored windows and sort globally ──────────────────
  fmt::print(stderr, "Step 4: Expanding and sorting windows...\n");
  auto windows = ExpandAndSortWindows(annotations, cache);
  const usize num_windows = windows.size();

  // ── Step 5: Pre-compute test sample set (by annotation) ─────────────────
  fmt::print(stderr, "Step 5: Building test sample set...\n");
  const auto test_set = BuildTestSampleSet(annotations);

  // ── Step 6: Set up queues and open outputs ──────────────────────────────
  const auto send_queue = std::make_shared<InputQueue>(num_windows);
  const auto recv_queue = std::make_shared<OutputQueue>(num_windows);
  const moodycamel::ProducerToken producer_token(*send_queue);

  for (usize idx = 0; idx < num_windows; ++idx) {
    send_queue->enqueue(producer_token, idx);
  }

  lancet::hts::BgzfOstream output_bed;
  if (!output_bed.Open(output_path, lancet::hts::BgzfFormat::BED)) {
    fmt::print(stderr, "Error: cannot open output: {}\n", output_path.string());
    return 1;
  }
  output_bed << fmt::format("#chrom\tstart\tend\tsource\tname\tregion\tscale\tregion_length\tscore\n");

  // Open test TSV
  std::ofstream test_out(test_path);
  test_out << "chrom\tstart\tend\tsource\tname\tregion\tscale\tregion_length\tscore\n";

  // ── Step 7: Launch worker threads ───────────────────────────────────────
  fmt::print(stderr, "Step 7: Launching {} workers for {} windows...\n\n", num_threads, num_windows);
  std::vector<std::jthread> workers;
  workers.reserve(num_threads);
  for (usize t = 0; t < num_threads; ++t) {
    workers.emplace_back(ScoringWorker, std::ref(producer_token),
                         send_queue, recv_queue, &windows, &annotations, &cache);
  }

  // ── Step 8: Main loop — collect results, flush in genome order ──────────
  std::vector<double> scores(num_windows, 0.0);
  absl::FixedArray<bool> done(num_windows);
  done.fill(false);

  usize num_completed = 0;
  usize idx_to_flush = 0;
  ScoredWindow incoming{};
  moodycamel::ConsumerToken consumer_token(*recv_queue);
  EtaTimer eta(num_windows);
  usize rows_written = 0;
  usize test_rows_written = 0;

  while (num_completed < num_windows) {
    if (!recv_queue->try_dequeue(consumer_token, incoming)) {
      using namespace std::chrono_literals;
      std::this_thread::sleep_for(100ms);
      continue;
    }

    scores[incoming.window_idx] = incoming.score;
    done[incoming.window_idx] = true;
    num_completed++;
    eta.Increment();

    if (num_completed % 500000 == 0 || num_completed == num_windows) {
      const auto elapsed = absl::FormatDuration(absl::Trunc(eta.Elapsed(), absl::Seconds(1)));
      const auto remaining = absl::FormatDuration(absl::Trunc(eta.EstimatedEta(), absl::Seconds(1)));
      fmt::print(stderr, "Progress: {:>7.2f}% | {}/{} | Elapsed: {} | ETA: {} @ {:.0f}/s\n",
                 eta.PercentDone(), num_completed, num_windows, elapsed, remaining,
                 eta.RatePerSecond());
    }

    // Flush completed windows in sorted order
    while (idx_to_flush < num_windows && done[idx_to_flush]) {
      const auto& win = windows[idx_to_flush];
      const auto& ann = annotations[win.parent_idx];
      const auto line = FormatLine(ann, win, scores[idx_to_flush]);

      output_bed << line;
      rows_written++;

      if (test_set.contains(win.parent_idx)) {
        test_out << line;
        test_rows_written++;
      }

      idx_to_flush++;
    }
  }

  // ── Step 9: Cleanup ─────────────────────────────────────────────────────
  std::ranges::for_each(workers, std::mem_fn(&std::jthread::request_stop));
  std::ranges::for_each(workers, std::mem_fn(&std::jthread::join));
  workers.clear();

  output_bed.Close();
  test_out.close();

  fmt::print(stderr,
             "\n=== Done ===\n"
             "Full output: {} ({} rows)\n"
             "Test output: {} ({} rows from {} sampled annotations)\n",
             output_path.string(), rows_written,
             test_path.string(), test_rows_written, test_set.size());

  return 0;
}
