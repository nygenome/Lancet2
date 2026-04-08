#include "lancet/core/window_builder.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <ios>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/cleanup/cleanup.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/match.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "absl/types/span.h"
#include "lancet/base/logging.h"
#include "lancet/base/types.h"
#include "lancet/core/window.h"
#include "spdlog/fmt/bundled/core.h"

extern "C" {
#include "htslib/hts.h"
#include "htslib/kstring.h"
}

namespace lancet::core {

WindowBuilder::WindowBuilder(const std::filesystem::path &ref_path, const Params &params)
    : mParams(params), mRefPtr(std::make_unique<hts::Reference>(ref_path)) {
  static constexpr usize DEFAULT_NUM_REGIONS_TO_ALLOCATE = 1024;
  mInputRegions.reserve(DEFAULT_NUM_REGIONS_TO_ALLOCATE);
}

void WindowBuilder::AddAllReferenceRegions() {
  static const auto should_exclude_chrom = [](std::string_view chrom) -> bool {
    return chrom == "MT" || chrom == "chrM" || absl::StartsWith(chrom, "GL") || absl::StartsWith(chrom, "chrUn") ||
           absl::StartsWith(chrom, "chrEBV") || absl::StartsWith(chrom, "HLA-") || absl::EndsWith(chrom, "_random") ||
           absl::EndsWith(chrom, "_alt") || absl::EndsWith(chrom, "_decoy");
  };

  const auto ref_chroms = mRefPtr->ListChroms();
  mInputRegions.reserve(mInputRegions.size() + ref_chroms.size());

  std::ranges::for_each(ref_chroms, [this](const hts::Reference::Chrom &chrom) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (should_exclude_chrom(chrom.Name())) return;
    this->mInputRegions.emplace_back(ParseRegionResult{.mChromName = chrom.Name(), .mRegionSpan = {1, chrom.Length()}});
  });
}

void WindowBuilder::AddRegion(const std::string &region_spec) {
  mInputRegions.emplace_back(mRefPtr->ParseRegion(region_spec.c_str()));
}

void WindowBuilder::AddBatchRegions(absl::Span<const std::string> region_specs) {
  mInputRegions.reserve(mInputRegions.size() + region_specs.size());
  std::ranges::for_each(region_specs, [this](const std::string &spec) {
    this->mInputRegions.emplace_back(this->mRefPtr->ParseRegion(spec.c_str()));
  });
}

void WindowBuilder::AddBatchRegions(const std::filesystem::path &bed_file) {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (bed_file.empty()) return;

  // Open the BED file using HTSlib's htsFile API instead of standard C++ streams.
  // This natively leverages the libcurl plugins built into Lancet2, allowing
  // direct streaming of BED files over s3://, gs://, and http(s):// protocols.
  htsFile *fp = hts_open(bed_file.c_str(), "r");
  if (fp == nullptr) {
    throw std::runtime_error(fmt::format("Could not open bed file: {}", bed_file.string()));
  }

  const absl::Cleanup stream_cleaner = [fp, &bed_file] {
    if (hts_close(fp) < 0) {
      LOG_WARN("Failed to properly close BED file stream: {}", bed_file.string());
    }
  };

  usize line_num = 0;
  i64 region_start = 0;
  i64 region_end = 0;

  std::string curr_chrom;
  std::vector<std::string_view> tokens;
  tokens.reserve(3);
  
  kstring_t line = KS_INITIALIZE;
  const absl::Cleanup kstring_cleaner = [&line] { ks_free(&line); };

  // Stream exactly line-by-line avoiding full-file buffering to minimize memory footprint
  while (hts_getline(fp, '\n', &line) > 0) {
    line_num++;
    const std::string_view line_view(line.s, line.l);

    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (line_view.starts_with('#') || line_view.empty()) continue;

    tokens = absl::StrSplit(line_view, absl::ByChar('\t'));

    if (tokens.size() != 3) {
      const auto msg = fmt::format("Invalid bed line with {} columns at line number {}", tokens.size(), line_num);
      throw std::runtime_error(msg);
    }

    if (!absl::SimpleAtoi(tokens[1], &region_start) || !absl::SimpleAtoi(tokens[2], &region_end)) {
      const auto msg = fmt::format("Could not parse line {} in bed: {}", line_num, bed_file.filename().string());
      throw std::runtime_error(msg);
    }

    curr_chrom = tokens[0];
    if (!mRefPtr->FindChromByName(curr_chrom).ok()) {
      const auto msg = fmt::format("Could not find chrom {} from bed file line {} in reference", tokens[0], line_num);
      throw std::runtime_error(msg);
    }

    mInputRegions.emplace_back(ParseRegionResult{.mChromName = curr_chrom, .mRegionSpan = {region_start, region_end}});
  }
}

auto WindowBuilder::StepSize(const Params &params) -> i64 {
  const auto val = (static_cast<f64>(100 - params.mPercentOverlap) / 100.0) * static_cast<f64>(params.mWindowLength);
  // round to ensure that steps always move in multiples of 100
  return static_cast<i64>(std::ceil(val / 100.0) * 100.0);
}

// ---------------------------------------------------------------------------
// ExpectedTargetWindows: arithmetic estimate without allocations
// ---------------------------------------------------------------------------

auto WindowBuilder::ExpectedTargetWindows() const -> usize {
  if (mInputRegions.empty()) return 0;

  const auto step_size = StepSize(mParams);
  const auto window_len = static_cast<i64>(mParams.mWindowLength);
  usize total_expected = 0;

  for (auto region : mInputRegions) {
    PadInputRegion(region);
    const auto region_len = static_cast<i64>(region.Length());

    if (region_len <= window_len) {
      total_expected += 1;
    } else {
      // Number of full windows that fit, plus account for the sliding step
      total_expected += static_cast<usize>((region_len - window_len) / step_size) + 1;
    }
  }

  return total_expected;
}

// ---------------------------------------------------------------------------
// SortInputRegions: pre-sort for deterministic sequential batch emission
// ---------------------------------------------------------------------------

void WindowBuilder::SortInputRegions() {
  // Deduplicate first
  std::ranges::sort(mInputRegions, [this](const ParseRegionResult &lhs, const ParseRegionResult &rhs) -> bool {
    const auto lhs_chrom = mRefPtr->FindChromByName(lhs.mChromName);
    const auto rhs_chrom = mRefPtr->FindChromByName(rhs.mChromName);
    const auto lhs_idx = lhs_chrom.ok() ? lhs_chrom->Index() : -1;
    const auto rhs_idx = rhs_chrom.ok() ? rhs_chrom->Index() : -1;
    if (lhs_idx != rhs_idx) return lhs_idx < rhs_idx;
    const auto lhs_start = lhs.mRegionSpan[0].value_or(0);
    const auto rhs_start = rhs.mRegionSpan[0].value_or(0);
    if (lhs_start != rhs_start) return lhs_start < rhs_start;
    return lhs.mRegionSpan[1].value_or(0) < rhs.mRegionSpan[1].value_or(0);
  });

  // Remove exact duplicates
  const auto last = std::unique(mInputRegions.begin(), mInputRegions.end());
  mInputRegions.erase(last, mInputRegions.end());
}

// ---------------------------------------------------------------------------
// BuildWindows: monolithic generation (for small region sets / targeted panels)
// ---------------------------------------------------------------------------

auto WindowBuilder::BuildWindows() const -> std::vector<WindowPtr> {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mInputRegions.empty()) return {};

  std::string rspec;
  const auto nregs = mInputRegions.size();
  const auto window_len = static_cast<i64>(mParams.mWindowLength);
  const auto pct_olap = static_cast<i64>(mParams.mPercentOverlap);
  LOG_INFO("Using {} input region(s) to build {}bp moving windows with {}% overlap", nregs, window_len, pct_olap)

  const auto step_size = StepSize(mParams);
  absl::flat_hash_set<WindowPtr> uniq_windows;

  for (ParseRegionResult region : mInputRegions) {
    PadInputRegion(region);
    const auto chrom = mRefPtr->FindChromByName(region.mChromName).value();

    if (region.Length() <= window_len) {
      auto wptr = std::make_shared<Window>(std::move(region), chrom, mRefPtr->FastaPath());
      uniq_windows.emplace(std::move(wptr));
      continue;
    }

    const auto chrom_has_colon = region.mChromName.find(':') != std::string::npos;
    auto curr_window_start = static_cast<i64>(region.mRegionSpan[0].value());
    const auto max_window_pos = static_cast<i64>(region.mRegionSpan[1].value());

    while ((curr_window_start + window_len) <= max_window_pos) {
      const i64 curr_window_end = curr_window_start + window_len;
      rspec = chrom_has_colon ? fmt::format("{{{}}}:{}-{}", region.mChromName, curr_window_start, curr_window_end)
                              : fmt::format("{}:{}-{}", region.mChromName, curr_window_start, curr_window_end);

      auto wptr = std::make_shared<Window>(mRefPtr->ParseRegion(rspec.c_str()), chrom, mRefPtr->FastaPath());
      uniq_windows.emplace(std::move(wptr));
      curr_window_start += step_size;
    }
  }

  static const auto WindowPtrComparator = [](const WindowPtr &first, const WindowPtr &second) -> bool {
    // NOLINTBEGIN(readability-braces-around-statements)
    if (first->ChromIndex() != second->ChromIndex()) return first->ChromIndex() < second->ChromIndex();
    if (first->StartPos1() != second->StartPos1()) return first->StartPos1() < second->StartPos1();
    return first->EndPos1() < second->EndPos1();
    // NOLINTEND(readability-braces-around-statements)
  };

  usize current_idx = 0;
  std::vector<WindowPtr> results(uniq_windows.cbegin(), uniq_windows.cend());
  std::ranges::sort(results, WindowPtrComparator);
  std::ranges::for_each(results, [&current_idx](WindowPtr &window_ptr) {
    window_ptr->SetGenomeIndex(current_idx);
    current_idx++;
  });

  return results;
}

// ---------------------------------------------------------------------------
// BuildWindowsBatch: pipelined emission for WGS memory control
//
// Requires SortInputRegions() to have been called first.
// Uses explicit reference parameters (region_idx, window_start, global_idx)
// to accurately track iteration state across repeated invocations, completely
// dropping the older complex offset encoding scheme.
// ---------------------------------------------------------------------------

auto WindowBuilder::BuildWindowsBatch(usize& region_idx, i64& window_start, usize& global_idx) const -> std::vector<WindowPtr> {
  if (mInputRegions.empty() || region_idx >= mInputRegions.size()) return {};

  const auto window_len = static_cast<i64>(mParams.mWindowLength);
  const auto step_size = StepSize(mParams);
  std::string rspec;

  std::vector<WindowPtr> batch;
  batch.reserve(BATCH_SIZE);

  while (region_idx < mInputRegions.size() && batch.size() < BATCH_SIZE) {
    auto region = mInputRegions[region_idx];
    PadInputRegion(region);
    
    const auto chrom = mRefPtr->FindChromByName(region.mChromName);
    if (!chrom.ok()) {
      region_idx++;
      window_start = -1;
      continue;
    }

    if (region.Length() <= window_len) {
      if (window_start == -1) {
        auto wptr = std::make_shared<Window>(std::move(region), chrom.value(), mRefPtr->FastaPath());
        wptr->SetGenomeIndex(global_idx++);
        batch.emplace_back(std::move(wptr));
        window_start = 1; // mark as done for this iteration
      }
      region_idx++;
      window_start = -1;
      continue;
    }

    const auto chrom_has_colon = region.mChromName.find(':') != std::string::npos;
    if (window_start == -1) {
      window_start = static_cast<i64>(region.mRegionSpan[0].value());
    }
    
    const auto max_window_pos = static_cast<i64>(region.mRegionSpan[1].value());

    while ((window_start + window_len) <= max_window_pos && batch.size() < BATCH_SIZE) {
      const i64 curr_window_end = window_start + window_len;
      rspec = chrom_has_colon
                  ? fmt::format("{{{}}}:{}-{}", region.mChromName, window_start, curr_window_end)
                  : fmt::format("{}:{}-{}", region.mChromName, window_start, curr_window_end);

      auto wptr = std::make_shared<Window>(mRefPtr->ParseRegion(rspec.c_str()), chrom.value(), mRefPtr->FastaPath());
      wptr->SetGenomeIndex(global_idx++);
      batch.emplace_back(std::move(wptr));
      
      window_start += step_size;
    }

    if ((window_start + window_len) > max_window_pos) {
      region_idx++;
      window_start = -1;
    }
  }

  return batch;
}

// ---------------------------------------------------------------------------
// PadInputRegion
// ---------------------------------------------------------------------------

void WindowBuilder::PadInputRegion(ParseRegionResult &result) const {
  const auto contig_info = mRefPtr->FindChromByName(result.mChromName);
  if (!contig_info.ok()) {
    const auto msg = fmt::format("No chromosome named {} found in reference", result.mChromName);
    throw std::runtime_error(msg);
  }

  const auto contig_max_len = contig_info.value().Length();
  const auto curr_start = result.mRegionSpan[0].value_or(1);
  const auto curr_end = result.mRegionSpan[1].value_or(contig_max_len);

  const auto start_underflows = curr_start <= mParams.mRegionPadding;
  const auto end_overflows = (curr_end > contig_max_len) || ((contig_max_len - curr_end) <= mParams.mRegionPadding);

  start_underflows ? result.mRegionSpan[0] = 1 : result.mRegionSpan[0] = (curr_start - mParams.mRegionPadding);
  end_overflows ? result.mRegionSpan[1] = contig_max_len : result.mRegionSpan[1] = (curr_end + mParams.mRegionPadding);

  if (result.Length() < mParams.mWindowLength) {
    const u64 diff = std::abs(static_cast<i64>(result.Length()) - static_cast<i64>(mParams.mWindowLength) - 1);
    const auto curr_left = result.mRegionSpan[0].value();
    const auto curr_right = result.mRegionSpan[1].value();
    const auto left_new_val = (diff / 2) > curr_left ? curr_left - 1 : curr_left - (diff / 2);
    const auto left_flank = curr_left - left_new_val;
    const auto goes_overmax = curr_right + (diff - left_flank) > contig_max_len;
    result.mRegionSpan[0] = curr_left - left_flank;
    result.mRegionSpan[1] = goes_overmax ? contig_max_len : curr_right + (diff - left_flank);
  }
}

}  // namespace lancet::core
